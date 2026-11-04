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
# (4) Figure 3 (sex differences in richness, diversity and composition)
# (5) Figure 4 (sex-biased genes in the gonad)
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 03_data_analysis/14_assemble_figures.R
#
# Reads:  data/03_data_analysis/05_overview_plots/plots/fig2[BCD]_*.rds
#         figures/schematics/fig2A_JS.png
#         figures/schematics/fig2E.png  (exported from manuscript/figures.pptx, slide 5)
#         data/03_data_analysis/05_overview_plots/plots/figS4_*.rds
#         data/03_data_analysis/12_isoformmetrics_bysex/plots/fig3*.rds, fig4*.rds
# Writes: figures/14_fig2.png, figures/14_fig2.pdf   (approved)
#         figures/14_figS4.png, figures/14_figS4.pdf (approved)
#         figures/14_fig3.png, figures/14_fig3.pdf   (approved)
#         figures/14_fig4.png, figures/14_fig4.pdf   (approved)
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

# macOS Quartz PDF keeps hyphens and fonts the same as the PNG (the default pdf device turns "-" into a long dash)
quartz_pdf <- function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...)

# Make the points of a saved plot bigger (multiplies every geom_point size)
bigger_points <- function(p, mult) {
  for (l in p$layers) {
    if (inherits(l$geom, "GeomPoint") && !is.null(l$aes_params$size)) {
      l$aes_params$size <- l$aes_params$size * mult
    }
  }
  p
}

# Make the text of a saved plot bigger (adds by_pt points to every geom_text size)
bigger_text <- function(p, by_pt) {
  for (l in p$layers) {
    if (inherits(l$geom, "GeomText") && !is.null(l$aes_params$size)) {
      l$aes_params$size <- l$aes_params$size + by_pt / .pt
    }
  }
  p
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

################################################################################
### (4) Figure 3
################################################################################

fig3_dir <- here::here("data", "03_data_analysis", "12_isoformmetrics_bysex", "plots")

female_col <- "#E87D72"
strip_size <- 11  # all Tail/Gonad strips the same size
male_col <- "#55BCC2"
# Text one point larger than the panel defaults (11 pt) everywhere in Fig 3
fig3_text <- theme(text = element_text(size = 12),
                   strip.text = element_text(size = strip_size + 1),
                   plot.tag = element_text(size = 21))

# Volcano (density on top): add a tissue strip and "more ... in females/males" arrows on the density
volcano_panel <- function(file, tissue, what) {
  v <- readRDS(file.path(fig3_dir, file))
  dens <- v[[1]]
  ymax <- max(layer_data(dens, 1)$y)
  y <- ymax * 0.55
  dens$data$strip <- tissue
  v[[2]] <- bigger_points(v[[2]], 1.3)
  v[[1]] <- dens +
    annotate("segment", x = -0.75, xend = -1.3, y = y, yend = y, colour = female_col,
             arrow = arrow(length = unit(0.08, "in")), linewidth = 0.6) +
    annotate("segment", x = 0.75, xend = 1.3, y = y, yend = y, colour = male_col,
             arrow = arrow(length = unit(0.08, "in")), linewidth = 0.6) +
    annotate("text", x = -1.4, y = y, label = paste0("more ", what, "\nin females"), size = 3.35, hjust = 1, lineheight = 0.9) +
    annotate("text", x = 1.4, y = y, label = paste0("more ", what, "\nin males"), size = 3.35, hjust = 0, lineheight = 0.9) +
    facet_wrap(~ strip) +
    theme(strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text = element_text(size = strip_size))
  v & fig3_text
}

# Arrow key beside the forest plots ("higher in males" up, "higher in females" down)
arrow_key <- function(up = "higher in\nmales", down = "higher in\nfemales",
                      up_col = male_col, down_col = female_col, text_size = 3.75) {
  ggplot() +
    annotate("segment", x = 0, xend = 0, y = 0.15, yend = 0.9, colour = up_col,
             arrow = arrow(length = unit(0.08, "in")), linewidth = 0.6) +
    annotate("segment", x = 0, xend = 0, y = -0.15, yend = -0.9, colour = down_col,
             arrow = arrow(length = unit(0.08, "in")), linewidth = 0.6) +
    annotate("text", x = 0.15, y = 0.5, label = up, hjust = 0, size = text_size, lineheight = 0.9) +
    annotate("text", x = 0.15, y = -0.5, label = down, hjust = 0, size = text_size, lineheight = 0.9) +
    coord_cartesian(xlim = c(-0.1, 1.2), ylim = c(-1, 1)) +
    theme_void()
}

tissue_labels <- labeller(tissue = c(tail = "Tail", gonad = "Gonad", Tail = "Tail", Gonad = "Gonad"))
forest_panel <- function(file) {
  readRDS(file.path(fig3_dir, file)) +
    facet_wrap(~ tissue, labeller = tissue_labels) +
    theme(strip.text = element_text(size = strip_size))
}

fig3A <- volcano_panel("fig3A_volcano_tail_richness.rds", "Tail", "isoforms") |
  volcano_panel("fig3A_volcano_gonad_richness.rds", "Gonad", "isoforms")
# Legend drawn below B's boxes without taking space, so rows A and D stay the same height
fig3B <- bigger_text(readRDS(file.path(fig3_dir, "fig3B_boxplot_richness.rds")), 1) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5, -0.2),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        strip.text = element_text(size = strip_size))
fig3C <- forest_panel("fig3C_composition.rds")
fig3D <- volcano_panel("fig3D_volcano_tail_diversity.rds", "Tail", "effective\nisoforms") |
  volcano_panel("fig3D_volcano_gonad_diversity.rds", "Gonad", "effective\nisoforms")
fig3E <- bigger_text(readRDS(file.path(fig3_dir, "fig3E_boxplot_diversity.rds")), 1) +
  theme(legend.position = "none", strip.text = element_text(size = strip_size))
fig3F <- forest_panel("fig3F_composition_effective.rds") +
  labs(y = "Change in effective\nisoform composition")

# One 2 x 4 grid so B/E and C/F line up between rows (volcanoes are wrapped, so they keep their own layout)
fig3 <- wrap_elements(full = fig3A) + fig3B + fig3C + arrow_key() +
  wrap_elements(full = fig3D) + fig3E + fig3F + arrow_key() +
  plot_layout(ncol = 4, widths = c(3.8, 0.6, 0.8, 0.4)) +  # wide volcanoes
  plot_annotation(tag_levels = list(c("A", "B", "C", "", "D", "E", "F", ""))) &
  tag_theme &
  fig3_text

# Approved 2026-11-04 (resized for readability): saved as the final figure
# Saved small so text and points stay readable when the figure is placed on the page
ggsave(here::here("figures", "14_fig3.png"),
       plot = fig3, width = 12.5, height = 6.6, dpi = 300)
ggsave(here::here("figures", "14_fig3.pdf"),
       plot = fig3, width = 12.5, height = 6.6, device = quartz_pdf)

################################################################################
### (5) Figure 4
################################################################################

# Uses fig3_dir, strip_size and arrow_key() from (4)
male_biased_col <- "#08519c"
female_biased_col <- "#E05A4E"

# Gonad strip on a panel without facets (tissue is gonad only)
gonad_strip <- function(p) {
  p$data$strip <- "Gonad"
  p + facet_wrap(~ strip) +
    theme(strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text = element_text(size = strip_size))
}
sex_labels <- labeller(sex = c(female = "Female", male = "Male"))
forest_sex_panel <- function(file) {
  readRDS(file.path(fig3_dir, file)) +
    facet_wrap(~ sex, labeller = sex_labels) +
    theme(strip.text = element_text(size = strip_size))
}

# A carries the shared legend (sex points + bias boxes), drawn between the rows without taking space
fig4A <- gonad_strip(readRDS(file.path(fig3_dir, "fig4A_boxplot_bias_richness.rds"))) +
  theme(legend.position = "inside",
        legend.position.inside = c(1.15, -0.17),
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank()) +
  # Sex points first, then bias boxes (the sex scale is the second colour scale)
  guides(colour = guide_legend(order = 2), colour_new = guide_legend(order = 1))
fig4B <- gonad_strip(readRDS(file.path(fig3_dir, "fig4B_fc_bias_richness.rds"))) +
  theme(legend.position = "none")
fig4C <- forest_sex_panel("fig4C_composition_bias.rds")
fig4D <- gonad_strip(readRDS(file.path(fig3_dir, "fig4D_boxplot_bias_diversity.rds"))) +
  theme(legend.position = "none")
fig4E <- gonad_strip(readRDS(file.path(fig3_dir, "fig4E_fc_bias_diversity.rds"))) +
  theme(legend.position = "none")
fig4F <- forest_sex_panel("fig4F_composition_bias_effective.rds")

bias_key <- arrow_key(up = "higher in\nmale-biased\ngenes", down = "higher in\nfemale-biased\ngenes",
                      up_col = male_biased_col, down_col = female_biased_col, text_size = 2.8)

fig4 <- fig4A + fig4B + fig4C + bias_key +
  fig4D + fig4E + fig4F + bias_key +
  plot_layout(ncol = 4, widths = c(1, 1, 1.25, 0.6)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "", "D", "E", "F", ""))) &
  tag_theme &
  theme(plot.margin = margin(2, 2, 2, 2))  # panels close together, like the preprint

# Approved 2026-11-03 ("will do for now"): saved as the final figure
ggsave(here::here("figures", "14_fig4.png"),
       plot = fig4, width = 8.5, height = 6.6, dpi = 300)
ggsave(here::here("figures", "14_fig4.pdf"),
       plot = fig4, width = 8.5, height = 6.6, device = quartz_pdf)
