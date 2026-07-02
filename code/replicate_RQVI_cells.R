library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# ── Paths ─────────────────────────────────────────────────────────────────────
data_path   <- "~/Desktop/Immgen/immgen-t-factors/data/cell_level_corr_rqvi_flashier.csv"
figure_path <- "~/Desktop/Immgen/immgen-t-factors/figures/Figure_RQVI/RQVI_CellLevel_Cor.pdf"

# ── Load & rename ─────────────────────────────────────────────────────────────
cell_level_corr_rqvi_flashier <- read.csv(data_path, header = TRUE) |>
  rename(flashier_GP = X) |>
  rename_with(~ paste0("RQVI_", seq_along(.)), .cols = -flashier_GP)

# ── Row-wise max correlation per GP ───────────────────────────────────────
rqvi_cols <- cell_level_corr_rqvi_flashier |> select(starts_with("RQVI_"))

cell_level_corr_rqvi_flashier <- cell_level_corr_rqvi_flashier |>
  mutate(
    max_corr = apply(rqvi_cols, 1, function(v) {
      m <- max((v), na.rm = TRUE)
      if (is.infinite(m)) NA_real_ else m
    }),
    well_replicated = max_corr > 0.5
  )

# ── Threshold sweep: how many GPs are well-replicated at each threshold? ──────
thr_df <- tibble(threshold = seq(0, 1, by = 0.01)) |>
  mutate(
    n_well    = vapply(threshold, function(t) {
      sum(cell_level_corr_rqvi_flashier$max_corr > t, na.rm = TRUE)
    }, integer(1)),
    prop_well = n_well / nrow(cell_level_corr_rqvi_flashier)
  )

ggplot(thr_df, aes(x = threshold, y = n_well)) +
  geom_line() +
  labs(x = "Threshold on max_corr", y = "# well-replicated GPs") +
  theme_minimal()

# ── Column-wise: which RQVI programs are not well replicated? ─────────────────
threshold <- 0.3

max_corr_perRQVI <- cell_level_corr_rqvi_flashier |>
  select(starts_with("RQVI_")) |>
  summarise(across(everything(), ~ max((.), na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "RQVI_program", values_to = "max_corr") |>
  mutate(well_replicated = max_corr > threshold)

sum(!max_corr_perRQVI$well_replicated)
max_corr_perRQVI$RQVI_program[!max_corr_perRQVI$well_replicated]

# ── Heatmap: GP rows vs RQVI columns ─────────────────────────────────────────
correlation_matrix <- read.csv(data_path, row.names = 1)
rownames(correlation_matrix) <- paste0("GP",   seq_len(nrow(correlation_matrix)))
colnames(correlation_matrix) <- paste0("RQVI_", seq_len(ncol(correlation_matrix)))

# Drop columns that are entirely NA
correlation_matrix <- correlation_matrix[
  , colSums(is.na(correlation_matrix)) < nrow(correlation_matrix)
]

pheatmap(
  mat            = correlation_matrix,
  show_rownames  = FALSE,
  show_colnames  = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  cluster_rows   = TRUE,
  color = c(
    colorRampPalette(c("#d0d8e4", "white"))(20),    # -0.8 to 0.0: pale blue-gray → white
    colorRampPalette(c("white", "#ffcccc"))(15),    #  0.0 to 0.3: white → very light pink
    colorRampPalette(c("#ffcccc", "#ff4444"))(15),  #  0.3 to 0.5: light pink → medium red
    colorRampPalette(c("#ff4444", "#cc0000"))(25),  #  0.5 to 0.7: medium → strong red
    colorRampPalette(c("#cc0000", "#4d0000"))(25)   #  0.7 to 1.0: strong → very dark red
  ),
  breaks = c(
    seq(-0.8, 0.0, length.out = 21),
    seq(0.0,  0.3, length.out = 16)[-1],
    seq(0.3,  0.5, length.out = 16)[-1],
    seq(0.5,  0.7, length.out = 26)[-1],
    seq(0.7,  1.0, length.out = 26)[-1]
  ),
  legend         = FALSE,
  filename       = figure_path,
  width          = 5,
  height         = 5
)
