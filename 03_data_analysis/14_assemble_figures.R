## Step 14: Assemble paper figures from saved panels
## Last modified by Linley Sherin 2026-11-03
#
# Runs locally. Builds the multi-panel paper figures with patchwork from the
# ggplot panels saved (.rds) by the analysis scripts, plus schematic images.
# New figures are written to figures/test/ until approved; approved ones go to figures/14_<fig>.png/.pdf.
#
# (1) Helpers
# (2) Figure 2
# (3) Supplementary Figure S4 (long- and short-read PCAs + heatmaps; reviewer)
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 03_data_analysis/14_assemble_figures.R
#
# Reads:  data/03_data_analysis/05_overview_plots/plots/fig2[BCD]_*.rds
#         figures/schematics/fig2A_JS.png
#         figures/schematics/fig2E.png  (exported from manuscript/figures.pptx, slide 5)
#         data/03_data_analysis/05_overview_plots/plots/figS4_*.rds
# Writes: figures/14_fig2.png, figures/14_fig2.pdf   (approved)
#         figures/14_figS4.png, figures/14_figS4.pdf (approved)
#
# Overwrites its outputs every time it runs.

library(ggplot2)
library(patchwork)
library(png)

################################################################################
### (1) Helpers
################################################################################

out_dir <- here::here("figures", "test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read a schematic image as a patchwork panel (keeps its aspect ratio).
# trim = TRUE cuts off blank (white or transparent) edges first.
image_panel <- function(path, trim = FALSE, pad_top = 0, pad_left = 0) {
  img <- png::readPNG(path)
  if (trim) {
    rgb_white <- apply(img[, , 1:3, drop = FALSE] > 0.98, c(1, 2), all)
    blank <- if (dim(img)[3] == 4) rgb_white | img[, , 4] < 0.02 else rgb_white
    rows <- which(!apply(blank, 1, all))
    cols <- which(!apply(blank, 2, all))
    img <- img[min(rows):max(rows), min(cols):max(cols), , drop = FALSE]
  }
  # Optional blank border (fraction of image height/width) so the panel letter sits outside the image
  if (pad_top > 0 || pad_left > 0) {
    nr <- dim(img)[1]; nc <- dim(img)[2]; nch <- dim(img)[3]
    top <- round(nr * pad_top); left <- round(nc * pad_left)
    padded <- array(1, dim = c(nr + top, nc + left, nch))
    padded[(top + 1):(top + nr), (left + 1):(left + nc), ] <- img
    img <- padded
  }
  wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE))
}

# Panel letters the same size and style in every figure
tag_theme <- theme(plot.tag = element_text(size = 20, face = "plain"))

################################################################################
### (2) Figure 2
################################################################################

fig2_dir <- here::here("data", "03_data_analysis", "05_overview_plots", "plots")

fig2A <- image_panel(here::here("figures", "schematics", "fig2A_JS.png"))
fig2B <- readRDS(file.path(fig2_dir, "fig2B_isoform_count_histogram.rds"))
fig2C <- readRDS(file.path(fig2_dir, "fig2C_relabund_histogram.rds"))
fig2D <- readRDS(file.path(fig2_dir, "fig2D_potential_codingstatus.rds")) +
  # Larger legend text; extra bottom margin makes the bars a bit shorter (closer to B and C)
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.8, "cm"),
        plot.margin = margin(t = 5.5, r = 5.5, b = 33, l = 5.5))
fig2E <- image_panel(here::here("figures", "schematics", "fig2E.png"), trim = TRUE, pad_top = 0.07, pad_left = 0.05)

# Top row: schematic A, histograms B and C. Bottom row: barplot D, schematic E
fig2_top <- fig2A + fig2B + fig2C + plot_layout(widths = c(4.94, 2.93, 2.93))
# D's bar area is a fixed 2.3 in so it keeps its shape; E takes the rest
fig2_bottom <- fig2D + fig2E + plot_layout(widths = unit(c(2.3, 1), c("in", "null")))

fig2 <- (fig2_top / fig2_bottom) +
  plot_layout(heights = c(0.8, 1)) +
  plot_annotation(tag_levels = "A") &
  tag_theme

# Approved 2026-11-03 (D shortened with a bottom margin): saved as the final figure
ggsave(here::here("figures", "14_fig2.png"),
       plot = fig2, width = 10.8, height = 6.3, dpi = 300)
# macOS Quartz PDF keeps hyphens and fonts the same as the PNG (the default pdf device turns "-" into a long dash)
quartz_pdf <- function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...)
ggsave(here::here("figures", "14_fig2.pdf"),
       plot = fig2, width = 10.8, height = 6.3, device = quartz_pdf)

################################################################################
### (3) Supplementary Figure S4
################################################################################

# PCAs (gene, isoform) and isoform heatmap for long and short reads
s4_dir <- here::here("data", "03_data_analysis", "05_overview_plots", "plots")

s4_pca_gene_long <- readRDS(file.path(s4_dir, "figS4_pca_gene_long.rds"))
s4_pca_iso_long <- readRDS(file.path(s4_dir, "figS4_pca_isoform_long.rds"))
s4_heat_long <- readRDS(file.path(s4_dir, "figS4_heatmap_long.rds")) +
  labs(title = "Isoform Expression (long reads)")
s4_pca_gene_short <- readRDS(file.path(s4_dir, "figS4_pca_gene_short.rds"))
s4_pca_iso_short <- readRDS(file.path(s4_dir, "figS4_pca_isoform_short.rds"))
s4_heat_short <- readRDS(file.path(s4_dir, "figS4_heatmap_short.rds")) +
  labs(title = "Isoform Expression (short reads)")

# Layout like the preprint (PCAs stacked beside the heatmap), long reads above short
s4_block <- function(pca_gene, pca_iso, heat) {
  (pca_gene / pca_iso | heat) + plot_layout(widths = c(1, 1.6))
}
figS4 <- (s4_block(s4_pca_gene_long, s4_pca_iso_long, s4_heat_long) /
              s4_block(s4_pca_gene_short, s4_pca_iso_short, s4_heat_short)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  tag_theme

# Approved 2026-11-03: saved as the final figure
ggsave(here::here("figures", "14_figS4.png"),
       plot = figS4, width = 11, height = 14, dpi = 300)
ggsave(here::here("figures", "14_figS4.pdf"),
       plot = figS4, width = 11, height = 14, device = quartz_pdf)
