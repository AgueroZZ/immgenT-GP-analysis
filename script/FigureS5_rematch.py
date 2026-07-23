"""Figure S5 rematch step. Re-derive the one-to-one EBMF<->RQVI matching on OUR
basis using Tianze's exact method, then aggregate the matched RQVI programs onto
the display clusters.

Matching basis: ALL common cells = our L_pm_filtered cells that also have RQVI
loadings, WITHOUT a non-thymocyte filter (aligns the two loading matrices
directly). Clusters = our annotation_level2 (108, including the thymocyte
cluster). This is the user-preferred alignment; vs a non-thymocyte-only matching
it changes only ~5/200 pairs and leaves the display-cluster correlation identical.

Display basis: non-thymocyte cells only, our annotation_level2 (107 clusters) --
this is what Figure S5 shows.

Method (identical to Tianze's fig_ebmf_rqvi_level2_comparison.py --recompute-matches):
  1. z-score every factor's mean-loading profile across clusters,
  2. signed Pearson r = EBMF_z^T @ RQVI_z / n_clusters,
  3. drop constant-profile RQVI candidates,
  4. maximum-weight one-to-one assignment via scipy linear_sum_assignment.
We only re-apply this method to the re-aligned data; the algorithm is unchanged.

EBMF loadings are read from Tianze's ebmf_cell_loadings.h5ad, which was verified
cell-for-cell identical to our L_pm_filtered (F_k == GP_k, cell-level r = 1.0).

Run: <miniforge>/envs/pyenv/bin/python script/FigureS5_rematch.py
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import linear_sum_assignment

FIG_DIR = Path("figures/generated/Figure S5")
PKG = Path("data/rqvi_loading/RQVI_EBMF_heatmap_data_v1/data")
EBMF_H5AD = PKG / "ebmf_cell_loadings.h5ad"
ALL_H5AD = PKG / "rqvi_all_10seeds_cell_loadings.h5ad"
# cell -> annotation table written by FigureS5.R (run step 1 first)
CELL_META = FIG_DIR / "S5_cell_metadata.csv.gz"
FACTORS = [f"F{k}" for k in range(1, 201)]


def _h5_index(f, group):
    g = f[group]
    key = g.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return np.array([v.decode() if isinstance(v, (bytes, bytearray)) else str(v) for v in g[key][:]])


def _zscore(M):
    m = M - M.mean(0, keepdims=True)
    s = M.std(0, ddof=0, keepdims=True)
    info = s.ravel() > np.finfo(float).eps
    z = np.zeros_like(M)
    z[:, info] = m[:, info] / s[:, info]
    return z, info


def _cluster_means(mat, codes, n_clusters):
    """mat: (n_cells x n_features) dense ndarray or scipy sparse; returns (K x n_features)."""
    counts = np.bincount(codes, minlength=n_clusters).astype(np.float64)
    ind = sp.csr_matrix((np.ones(codes.size), (codes, np.arange(codes.size))),
                        shape=(n_clusters, codes.size))
    prod = ind @ mat
    if sp.issparse(prod):
        prod = np.asarray(prod.todense())
    return np.asarray(prod, dtype=np.float64) / counts[:, None], counts


def main() -> None:
    # our L_pm_filtered cells -> annotation (no filtering yet)
    meta = pd.read_csv(CELL_META, dtype=str)
    cell_l1 = dict(zip(meta["cellID"], meta["annotation_level1"]))
    cell_l2 = dict(zip(meta["cellID"], meta["annotation_level2"]))

    # EBMF (F1..F200, dense) and RQVI all-2560 (sparse), same cell order
    fe = h5py.File(EBMF_H5AD, "r")
    obs = _h5_index(fe, "obs")
    e_var = _h5_index(fe, "var")
    e_pos = {n: i for i, n in enumerate(e_var)}
    E = fe["X"][:][:, [e_pos[c] for c in FACTORS]]              # n_obs x 200
    fe.close()

    fr = h5py.File(ALL_H5AD, "r")
    r_obs = _h5_index(fr, "obs")
    r_var = _h5_index(fr, "var")
    g = fr["X"]
    R = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g.attrs["shape"]))
    fr.close()
    if not np.array_equal(obs, r_obs):
        raise ValueError("EBMF and RQVI cell orders differ")

    # ---- MATCHING basis: all common cells (no non-thymocyte filter) ----
    m_mask = np.array([c in cell_l2 for c in obs])
    m_idx = np.where(m_mask)[0]
    m_labs = np.array([cell_l2[obs[i]] for i in m_idx])
    m_clusters = sorted(set(m_labs))
    m_code = {l: i for i, l in enumerate(m_clusters)}
    m_codes = np.array([m_code[l] for l in m_labs])
    K_match = len(m_clusters)

    ebmf_m, _ = _cluster_means(E[m_idx], m_codes, K_match)      # K_match x 200
    rqvi_m, _ = _cluster_means(R[m_idx], m_codes, K_match)      # K_match x 2560
    ez, e_info = _zscore(ebmf_m)
    rz, r_info = _zscore(rqvi_m)
    if not e_info.all():
        raise ValueError("constant EBMF factor on matching basis")
    corr = ez.T @ rz / K_match
    corr[:, ~r_info] = -np.inf
    cost = np.where(np.isfinite(corr), corr, -1e9)
    rows, cols = linear_sum_assignment(-cost)
    if rows.size != 200 or np.unique(cols).size != 200:
        raise RuntimeError("assignment did not cover all 200 EBMF factors uniquely")
    sel = np.empty(200, dtype=int)
    sel[rows] = cols
    match_r = corr[np.arange(200), sel]
    matched_candidates = r_var[sel]
    print(f"MATCHING basis: {len(m_idx)} common cells, {K_match} level2 clusters "
          f"(incl. {sorted(set(m_clusters) - set(cell_l2[c] for c in obs if c in cell_l2 and cell_l1[c] != 'thymocyte'))})")

    # ---- DISPLAY basis: non-thymocyte cells, order from S5_cluster_order.csv ----
    order = pd.read_csv(FIG_DIR / "S5_cluster_order.csv").sort_values("display_column")
    disp_clusters = order["level2_cluster"].astype(str).tolist()
    d_code = {l: i for i, l in enumerate(disp_clusters)}
    d_mask = np.array([(c in cell_l2) and (cell_l1[c] != "thymocyte") for c in obs])
    d_idx = np.where(d_mask)[0]
    d_codes = np.array([d_code[cell_l2[obs[i]]] for i in d_idx])
    rqvi_disp_all, d_counts = _cluster_means(R[d_idx], d_codes, len(disp_clusters))   # 107 x 2560
    if not np.array_equal(d_counts.astype(int), order["n_cells"].to_numpy()):
        raise ValueError("display cell counts differ from S5_cluster_order")
    rematched_disp = rqvi_disp_all[:, sel]                        # 107 x 200

    out = pd.DataFrame(rematched_disp, index=disp_clusters, columns=FACTORS)
    out.index.name = "level2_cluster"
    out.to_csv(FIG_DIR / "S5_rqvi_rematched_raw_means_level2.csv")

    # display-basis correlation for the caption
    ebmf_disp, _ = _cluster_means(E[d_idx], d_codes, len(disp_clusters))
    ez_d = _zscore(ebmf_disp)[0]
    rz_d = _zscore(rematched_disp)[0]
    disp_r = (ez_d * rz_d).sum(0) / len(disp_clusters)

    shipped = pd.read_csv(PKG / "ebmf_rqvi_multiseed_level2_one_to_one_matches.csv")
    shipped_map = dict(zip(shipped["ebmf_factor"].astype(str), shipped["rqvi_candidate"].astype(str)))
    same = sum(matched_candidates[k] == shipped_map.get(f"F{k+1}") for k in range(200))

    pd.DataFrame({
        "ebmf_factor": FACTORS,
        "rematched_rqvi_candidate": matched_candidates,
        "pearson_r_match_basis_108": match_r,
        "pearson_r_display_basis_107": disp_r,
        "same_as_shipped": [matched_candidates[k] == shipped_map.get(f"F{k+1}") for k in range(200)],
    }).to_csv(FIG_DIR / "S5_rematch_mapping.csv", index=False)

    print("=== rematch (matching basis = all common cells, 108 clusters) ===")
    print(f"match-basis r:   median {np.median(match_r):.3f}, >=0.5 {100*np.mean(match_r>=0.5):.1f}%, "
          f"r<0.3 {int((match_r<0.3).sum())}")
    print(f"display-basis r: median {np.median(disp_r):.3f}, >=0.5 {100*np.mean(disp_r>=0.5):.1f}%, "
          f"r<0.3 {int((disp_r<0.3).sum())}")
    print(f"pairs identical to Tianze's shipped matching: {same}/200")
    print(f"wrote {FIG_DIR/'S5_rqvi_rematched_raw_means_level2.csv'} and S5_rematch_mapping.csv")


if __name__ == "__main__":
    main()
