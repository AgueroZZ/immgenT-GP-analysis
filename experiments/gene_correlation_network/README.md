# Gene–gene correlation network

**Goal:** visualize relationships *between genes* to motivate gene programs (GPs)
— show that when you correlate each gene's GP-loading profile and connect
strongly-correlated genes, the emergent network structure recapitulates the GPs
(e.g. Treg genes around `Foxp3` form their own community).

All inputs come from the flashier factor matrix `data/F_pm_filtered.rds`
(19805 genes × 200 GPs), always **per-GP (per-column) max-abs normalized**.

## → Final figure

**`gp_gene_bipartite_up.png` / `.pdf`** — Fig-3d-style bipartite over ~all 200 GPs,
**up-regulated genes only** (default top-10, |loading|≥0.1):

- GPs = uniform blue balls (numbered); shared genes = small tan balls.
- Each marker gene (default `Foxp3, Cd8a, Cd8b1, Izumo1r, Il2ra`) has its own
  highlight color shared by its node and its (thicker) edges, listed in a side
  **color legend** — no on-plot gene text, so no label overlap; each marker's GP
  "territory" is traceable by edge color.
- Background (non-marker) edges are thin dark lines so the network structure
  reads without drowning the colored marker edges.
- FR layout with **radial compression** (r^0.65) to pull the sparse periphery in
  and reclaim whitespace; coords jittered slightly to separate structural-twin
  GPs; all GP numbers shown. Private (single-GP) gene dots hidden; single-GP
  star components dropped.

Built by **`gp_gene_bipartite_up.R [n_top] [drop_df]`** — edit the `MARKERS`
vector to change which genes are highlighted. `drop_df` = optional filter
dropping "glue" genes shared by more than `drop_df` GPs (markers exempt);
**default `Inf` = off** (pure top-`n_top` ≥ 0.1 rule). This is the version wired
into `script/Figure1.R` as panel 1C.

## → GP-highlight variant

**`gp_highlight_formal.png` / `.pdf`** (and `gp_highlight_internal.*`) — same
GP-gene bipartite network, but highlighting a chosen set of **GPs** instead of
marker genes. Each highlighted GP's color flows to its node and its edges, and
the genes it connects to are labeled; a chosen GP's signature is read off
directly and overlaps between GPs are visible.

- **formal** = no legend, no GP-index labels (for figures — the legend squeezes
  the plot).
- **internal** = legend + GP-index labels (for the team, to see color → GP).
- top-5 up genes per GP, plain FR layout (seed 1), large canvas for readable
  labels.

Built by **`gp_highlight_selected.R`** — refactored into `build_gp_gene_graph()`
/ `add_gp_highlights()` / `make_gp_plot()`; **to highlight a different GP set,
edit the `GP_HIGHLIGHTS` (GP → color) vector at the top and rerun. Nothing else
needs to change.**

## → Loading heatmap (now Figure 1D)

**`heatmap_loading.png` / `.pdf`** — giant heatmap of all 200 GP loadings across a
stratified cell sample (rows = cells ordered by lineage × organ, columns = GPs
clustered by similarity). The **same GPs highlighted in the GP-gene network** are
marked in the matching colors — a top annotation bar, colored/bold column labels,
and a box around each highlighted GP **column** (GPs are columns here). Built by
**`heatmap_loading.R`** (edit the `GP_HIGHLIGHTS` map to match the network).

This is the **internal version** (all 200 GP indices labeled + a color→GP
legend), for collaborators. The **publication version** used in the paper as
**Figure 1D** (`script/Figure1.R`, panel 1D) is the same heatmap but with the
background GP indices dropped and the highlighted-GP legend removed — only the
highlighted GPs are labeled.

Everything else is the exploration trail in **`archive/`** — kept for the record
of what was tried and why (see Findings below), not needed to read the result.

## Exploration trail (`archive/`)

| script | what it does |
|---|---|
| `explore_scale.R` | Gauge edge counts / degree at corr thresholds over **all 19805 genes**. |
| `build_network.R` | Full 19805-gene graph at corr ≥ 0.5, DrL layout. → hairball. |
| `signature_communities.R` | Signature genes (n=8053), corr ≥ 0.5, Louvain + Foxp3/Treg diagnostic. |
| `signature_layout.R` | FR layout of the corr ≥ 0.5 signature graph. → collapses. |
| `knn_network.R [k] [floor]` | kNN backbone, Louvain + FR layout. |
| `umap_network.R` | Gene UMAP + kNN edges. Clean but ≈ existing t-SNE atlas. |
| `gp_network.R [min_shared]` | GP-centric v1: GP nodes, edge = raw shared-gene count. |
| `gp_network_specific.R [df_max] [min_w]` | GP-centric v2: specificity-weighted edges. |
| `gp_gene_bipartite_all.R [drop_df]` | Fig-3d-style bipartite for all 200 GPs, up+down genes. Superseded by the up-only final figure. |
| `gp_network_final.R [df_max] [min_w]` | GP-GP network with Louvain-family convex hulls, node size ∝ PVE. Alternative broad-picture form (not the chosen final). |

## Findings so far

1. **corr ≥ 0.5 over all genes is a near-complete hairball**: 1.37M edges,
   99.7% of genes connected, median degree 108. Broadly-expressed housekeeping
   genes (`Nipbl`, `Cnot1`, `Psmd7`, `Cct7`…) correlate with almost everything
   and glue all clusters together. DrL/FR force layouts collapse it into one
   blob — no separable structure. (`build_network.R`, `signature_layout.R`)

2. **Community structure IS real** among signature genes: Louvain modularity
   0.82, and several communities are near-pure single GPs (GP75 92%, GP175 97%,
   GP92 94%). (`signature_communities.R`)

3. **The Treg program is a small, sparse module.** `Foxp3` has only ~13 strong
   correlates, several of them Treg genes (`Il2ra`, `Ikzf4`, `Nrp1`, `Lrrc32`,
   `Gpr83`). At default Louvain resolution it is absorbed into a giant
   housekeeping community; only at higher resolution (~8) does it resolve into a
   clean 16-gene, GP68-dominated Treg community. → the network is **multi-scale**.

4. **A global correlation threshold is the wrong edge rule for layout.** The
   kNN backbone (top-k correlates per gene) keeps the graph sparse and adaptive,
   so small specific programs keep their own edges instead of being glued to the
   housekeeping mass — force layouts then separate communities cleanly.

5. **Gene-level layout is a dead end; go GP-centric.** Force layouts (FR/DrL)
   collapse the small-world gene graph into a blob, and a UMAP of gene loadings
   just reproduces the existing t-SNE atlas. Switching nodes to **GPs** (≤200,
   layout-able, denoised) gives a genuinely structured community plot.

6. **Specificity weighting is what makes GP families resolve.** Linking GPs by
   raw shared-top-gene count glues all immune GPs together (broadly-expressed
   Ccl5 / Tmsb4x / granzymes / ribosomal sit in dozens of top-sets;
   modularity 0.39). Weighting each shared gene by 1/df and dropping df>10 glue
   genes lifts modularity to **0.82** and yields ~20 interpretable families:
   hepatocyte/intestinal (Alb, Apoa4, Fabp2), placental (Prl3b1, Tpbpa),
   pancreatic acinar (Ctrb1, Prss2), Paneth (Defa), stromal (Cryab, A2m, Dcn),
   IFN (Isg15, Stat1, Irf1), Th17/ILC (Il17a, Il1rl1, Il12rb2), cell cycle
   (Mcm3, Cenpa, histones), and a T-cell development family containing Foxp3/GP68
   (Satb1, Tox, Maf, Pdcd1, Themis, Ikzf2). → `gp_network_specific.R`

## Convention notes
- Dominant GP = `max.col` of the per-GP-normalized loading (original GP index),
  colored on the same 200-hue rainbow as the gene t-SNE atlas
  (`experiments/gene_umap_gp_space/`).
- All plots rasterized (`ggrastr`) for manageable PDF size.
