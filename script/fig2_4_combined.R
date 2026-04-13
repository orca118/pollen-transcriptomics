library(ggplot2)
library(cowplot)
library(png)
library(grid)

fig_dir <- file.path(getwd(), "figures")

read_img <- function(path) {
  png::readPNG(path)
}

make_panel <- function(path) {
  ggdraw() + draw_image(read_img(path))
}

# ============================
# FIGURE 2 (balanced panels)
# ============================
fig2 <- plot_grid(
  make_panel(file.path(fig_dir, "Fig2A_bulk_volcano.png")),
  make_panel(file.path(fig_dir, "Fig2B_bulk_GSEA.png")),
  labels = c("A", "B"),
  label_size = 18,
  ncol = 2,
  rel_widths = c(1, 1)   # equal size
)

ggsave(
  file.path(fig_dir, "Fig2_combined.png"),
  fig2,
  width = 10,   # 🔥 narrower = less empty space
  height = 5,
  dpi = 300,
  bg = "white"
)

# ============================
# FIGURE 3 (UMAP slightly larger)
# ============================
fig3 <- plot_grid(
  make_panel(file.path(fig_dir, "Fig3A_scRNA_UMAP_celltypes.png")),
  make_panel(file.path(fig_dir, "Fig3B_scRNA_cell_composition.png")),
  labels = c("A", "B"),
  label_size = 18,
  ncol = 2,
  rel_widths = c(1.3, 1)   # 🔥 UMAP slightly larger
)

ggsave(
  file.path(fig_dir, "Fig3_combined.png"),
  fig3,
  width = 11,
  height = 5,
  dpi = 300
)

# ============================
# FIGURE 4 (tight vertical layout)
# ============================
fig4 <- plot_grid(
  make_panel(file.path(fig_dir, "Fig4A_key_immune_genes_heatmap.png")),
  make_panel(file.path(fig_dir, "Fig4B_pathway_summary_heatmap.png")),
  labels = c("A", "B"),
  label_size = 18,
  ncol = 1,
  rel_heights = c(1.1, 1)   # 🔥 top slightly larger
)

ggsave(
  file.path(fig_dir, "Fig4_combined.png"),
  fig4,
  width = 8,    # 🔥 narrower = cleaner
  height = 12,
  dpi = 300
)

library(ggplot2)
library(cowplot)
library(png)
library(grid)

fig_dir <- file.path(getwd(), "figures")

read_img <- function(path) {
  png::readPNG(path)
}

make_panel <- function(path) {
  ggdraw() + draw_image(read_img(path))
}

# ============================
# FIGURE 4 (side-by-side)
# ============================
fig4_side <- plot_grid(
  make_panel(file.path(fig_dir, "Fig4A_key_immune_genes_heatmap.png")),
  make_panel(file.path(fig_dir, "Fig4B_pathway_summary_heatmap.png")),
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1.15, 1)   # left panel slightly wider
)

ggsave(
  file.path(fig_dir, "Fig4_combined_side_by_side.png"),
  fig4_side,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)