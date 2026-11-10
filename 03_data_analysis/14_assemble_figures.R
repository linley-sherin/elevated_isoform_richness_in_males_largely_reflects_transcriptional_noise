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
# (6) Figure 5 (isoform switching)
# (7) Figure 1 (long vs short reads; isoform- vs event-level)
#
# Usage (from the project folder):
#   /usr/local/bin/Rscript 03_data_analysis/14_assemble_figures.R
#
# Reads:  data/03_data_analysis/05_overview_plots/plots/fig2[BCD]_*.rds
#         figures/schematics/fig2A_JS.png
#         figures/schematics/fig2E.png  (exported from manuscript/figures.pptx, slide 5)
#         data/03_data_analysis/05_overview_plots/plots/figS4_*.rds
#         data/03_data_analysis/12_isoformmetrics_bysex/plots/fig3*.rds, fig4*.rds
#         data/03_data_analysis/13_isoformswitch/plots/fig5*.rds
#         figures/schematics/fig5A_JS.png
#         data/03_data_analysis/11_suppa2_parse/plots/fig1[BEF]_*.rds
#         data/03_data_analysis/07_data_correlation/plots/fig1C_corr.rds
# Writes: figures/14_fig2.png, figures/14_fig2.pdf   (approved)
#         figures/14_figS4.png, figures/14_figS4.pdf (approved)
#         figures/14_fig3.png, figures/14_fig3.pdf   (approved)
#         figures/14_fig4.png, figures/14_fig4.pdf   (approved)
#         figures/14_fig5.png, figures/14_fig5.pdf   (approved)
#         figures/test/14_fig1.png, figures/test/14_fig1.pdf (not yet approved)
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
  v[[2]] <- bigger_points(v[[2]], 1.3) + fig3_axis
  v[[1]] <- dens +
    annotate("segment", x = -0.75, xend = -1.3, y = y, yend = y, colour = female_col,
             arrow = arrow(length = unit(0.08, "in")), linewidth = 0.6) +
    annotate("segment", x = 0.75, xend = 1.3, y = y, yend = y, colour = male_col,
             arrow = arrow(length = unit(0.08, "in")), linewidth = 0.6) +
    # 11.2 pt = 6 pt at the journal's 170 mm width
    annotate("text", x = -1.4, y = y, label = "more\nin females", size = 11.2 / .pt, hjust = 1, lineheight = 0.9) +
    annotate("text", x = 1.4, y = y, label = "more\nin males", size = 11.2 / .pt, hjust = 0, lineheight = 0.9) +
    facet_wrap(~ strip) +
    theme(strip.background = element_rect(fill = "grey90", colour = NA),
          strip.text = element_text(size = strip_size),
          axis.text = element_blank(), axis.ticks = element_blank())  # density axes hidden
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

# Bigger tick labels (>= 6 pt at the journal's 170 mm width); added per panel so hidden axes stay hidden
fig3_axis <- theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

tissue_labels <- labeller(tissue = c(tail = "Tail", gonad = "Gonad", Tail = "Tail", Gonad = "Gonad"))
forest_panel <- function(file) {
  readRDS(file.path(fig3_dir, file)) +
    facet_wrap(~ tissue, labeller = tissue_labels) +
    theme(strip.text = element_text(size = strip_size), axis.text.y = element_text(size = 12))
}

fig3A <- volcano_panel("fig3A_volcano_tail_richness.rds", "Tail", "isoforms") |
  volcano_panel("fig3A_volcano_gonad_richness.rds", "Gonad", "isoforms")
# Legend drawn below B's boxes without taking space, so rows A and D stay the same height
fig3B <- bigger_text(readRDS(file.path(fig3_dir, "fig3B_boxplot_richness.rds")), 1) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5, -0.2),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        strip.text = element_text(size = strip_size),
        axis.text.y = element_text(size = 12))
fig3C <- forest_panel("fig3C_composition.rds")
fig3D <- volcano_panel("fig3D_volcano_tail_diversity.rds", "Tail", "effective\nisoforms") |
  volcano_panel("fig3D_volcano_gonad_diversity.rds", "Gonad", "effective\nisoforms")
fig3E <- bigger_text(readRDS(file.path(fig3_dir, "fig3E_boxplot_diversity.rds")), 1) +
  theme(legend.position = "none", strip.text = element_text(size = strip_size),
        axis.text.y = element_text(size = 12))
fig3F <- forest_panel("fig3F_composition_effective.rds") +
  labs(y = "Change in effective\nisoform composition")

# One 2 x 4 grid so B/E and C/F line up between rows (volcanoes are wrapped, so they keep their own layout)
fig3 <- wrap_elements(full = fig3A) + fig3B + fig3C + arrow_key() +
  wrap_elements(full = fig3D) + fig3E + fig3F + arrow_key() +
  plot_layout(ncol = 4, widths = c(3.8, 0.6, 0.8, 0.4)) +  # wide volcanoes
  plot_annotation(tag_levels = list(c("A", "B", "C", "", "D", "E", "F", ""))) &
  tag_theme &
  fig3_text

# Approved 2026-11-04 (resized for readability); text >= 6 pt at 170 mm approved 2026-11-10
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
          strip.text = element_text(size = strip_size),
          axis.text = element_blank(), axis.ticks = element_blank())  # density axes hidden
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

################################################################################
### (6) Figure 5
################################################################################

fig5_dir <- here::here("data", "03_data_analysis", "13_isoformswitch", "plots")
fig5_bg <- "#F8F4EF"  # beige box behind A
fig5_text <- theme(text = element_text(size = 13),
                   axis.text = element_text(size = 12))

# A: switching genes per switch type, rows in the schematic's order
# (male-biased count : female-biased count, top to bottom)
switch_rows <- c("1Female:1Male", "MultiFemale:1Male", "1Female:MultiMale", "MultiFemale:MultiMale")
fig5A_saved <- readRDS(file.path(fig5_dir, "fig5A_switching.rds"))
fig5A_data <- fig5A_saved$data
fig5A_data$row <- 5 - match(fig5A_data$switch_type, switch_rows)  # top row = 4
switch_cols <- fig5A_saved$scales$get_scales("fill")$palette(5)

# Schematic rows span 27.6-95.4% of the image height; scale it so they line up with bar rows 4-1
icons <- png::readPNG(here::here("figures", "schematics", "fig5A_JS.png"))
units_per_img <- 4 / (0.954 - 0.276)
img_top <- 4.5 + 0.276 * units_per_img
img_bottom <- img_top - units_per_img
icon_x <- -205  # left edge of the schematic in bar units

fig5A <- ggplot(fig5A_data, aes(x = Count, y = row, fill = orf_outcome)) +
  annotation_custom(grid::rasterGrob(icons, x = 0, hjust = 0, interpolate = TRUE),
                    xmin = icon_x, xmax = 0, ymin = img_bottom, ymax = img_top) +
  geom_col(width = 0.9, orientation = "y", position = position_stack(reverse = TRUE)) +
  annotate("text", x = 100, y = 4.95, label = "Isoform-switching genes", size = 14 / .pt) +
  scale_fill_manual(values = switch_cols, name = NULL) +
  scale_x_continuous(breaks = seq(0, 150, 50), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # crop the schematic's blank top and bottom; clip = "off" still draws it
  coord_cartesian(xlim = c(icon_x, 210), ylim = c(0.45, img_top - 0.25), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void() +
  fig5_text +
  theme(axis.text.x = element_text(size = 12, margin = margin(t = 4)),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.55, "cm"),
        plot.background = element_rect(fill = fig5_bg, colour = NA),
        plot.margin = margin(8, 8, 8, 8))

# B, C: unbiased isoforms and loss of coding (from script 13)
fig5_bar <- function(file) {
  p <- readRDS(file.path(fig5_dir, file)) + fig5_text + theme(plot.margin = margin(12, 6, 2, 6))
  # counts above the bars have no size set in 13 (default); make them 14 pt
  for (l in p$layers) {
    if (inherits(l$geom, "GeomText") && is.null(l$aes_params$size)) l$aes_params$size <- 14 / .pt
  }
  p
}
fig5B <- fig5_bar("fig5B_unbiased.rds") +
  scale_x_discrete(labels = function(x) sub(" \n", "\n", sub("unbiased ", "unbiased\n", x)))  # 3-line labels so they don't touch
fig5C <- fig5_bar("fig5C_loss.rds")

# wrap A so its legend does not squeeze B and C
fig5 <- wrap_elements(full = fig5A) / (fig5B | fig5C) +
  plot_layout(heights = c(1.3, 0.8)) +
  plot_annotation(tag_levels = "A") &
  tag_theme

# Approved 2026-11-10 (B/C row shortened, bar labels bigger): saved as the final figure
ggsave(here::here("figures", "14_fig5.png"),
       plot = fig5, width = 8.5, height = 7.4, dpi = 300)
ggsave(here::here("figures", "14_fig5.pdf"),
       plot = fig5, width = 8.5, height = 7.4, device = quartz_pdf)

################################################################################
### (7) Figure 1
################################################################################

fig1_dir <- here::here("data", "03_data_analysis", "11_suppa2_parse", "plots")
event_cols <- c(ES = "#94CCEC", MX = "#ECB4A0", A3A5 = "#9DBCB5", IR = "#8676BF")
exon_grey <- "grey75"
fig1_text <- theme(text = element_text(size = 12), axis.text = element_text(size = 11))

# Schematic helpers: exons as boxes on a line, splice junctions as peaked arcs
exons <- function(df, y, fill = exon_grey) {
  geom_rect(data = df, aes(xmin = start, xmax = end, ymin = y - 0.25, ymax = y + 0.25),
            fill = fill, inherit.aes = FALSE)
}
junction <- function(x1, x2, y, up = TRUE, h = 0.5) {
  ys <- if (up) y + 0.25 else y - 0.25
  yp <- if (up) ys + h else ys - h
  annotate("path", x = c(x1, (x1 + x2) / 2, x2), y = c(ys, yp, ys), colour = "grey35", linewidth = 0.35)
}
backbone <- function(x1, x2, y) annotate("segment", x = x1, xend = x2, y = y, yend = y, colour = "grey35", linewidth = 0.35)
schematic_theme <- theme_void() + theme(plot.title = element_text(size = 13), plot.margin = margin(5, 5, 5, 30))  # left margin keeps the title clear of the panel letter

# A: isoform-level (full-length transcripts of one gene)
iso <- data.frame(row = rep(4:1, c(4, 4, 3, 3)),
                  start = c(0, 5.2, 10.4, 15.6, 0, 5.2, 11.2, 15.6, 2.6, 11.2, 15.6, 2.6, 7.8, 11.2),
                  end   = c(3.6, 8.8, 14.0, 19.2, 3.6, 9.6, 14.0, 19.2, 6.2, 14.0, 19.2, 6.2, 9.7, 15.6))
iso_introns <- do.call(rbind, lapply(split(iso, iso$row), function(d)
  data.frame(row = d$row[-1], x1 = head(d$end, -1), x2 = d$start[-1])))
fig1A <- ggplot() +
  geom_segment(data = iso_introns, aes(x = x1, xend = x2, y = row, yend = row), colour = "grey60", linewidth = 0.8) +
  geom_rect(data = iso, aes(xmin = start, xmax = end, ymin = row - 0.17, ymax = row + 0.17), fill = "grey60") +
  labs(title = expression(underline("Isoform-level"))) +
  coord_cartesian(xlim = c(-1.5, 20.7), ylim = c(0, 4.5), expand = FALSE) +  # white space either side
  schematic_theme

# D: event-level (local splicing events), one schematic per event type
event_panel <- function(title, draw) {
  ggplot() + draw + labs(title = title) + coord_cartesian(xlim = c(0, 10), ylim = c(-1.1, 1.1)) +
    theme_void() + theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = 0)))
}
ex <- function(x1, x2, fill = exon_grey, y = 0) annotate("rect", xmin = x1, xmax = x2, ymin = y - 0.25, ymax = y + 0.25, fill = fill)
d_es <- event_panel("Exon skipping (ES)", list(backbone(1, 9, 0), ex(1.5, 3.5), ex(4.2, 5.8, event_cols["ES"]), ex(6.5, 8.5),
  junction(3.5, 4.2, 0), junction(5.8, 6.5, 0), junction(3.5, 6.5, 0, up = FALSE)))
d_mx <- event_panel("Mutually exclusive exons (MX)", list(backbone(0.5, 9.5, 0), ex(1, 2.6), ex(3.3, 4.7, event_cols["MX"]),
  ex(5.4, 6.8, event_cols["MX"]), ex(7.5, 9), junction(2.6, 3.3, 0), junction(4.7, 7.5, 0),
  junction(2.6, 5.4, 0, up = FALSE), junction(6.8, 7.5, 0, up = FALSE)))
d_a35 <- event_panel("Alternative 3\u2032/5\u2032 splice sites (A3/A5)", list(
  backbone(0, 4.6, 0), ex(0.4, 1.8), ex(1.8, 2.3, event_cols["A3A5"]), ex(3.2, 4.4),
  junction(2.3, 3.2, 0), junction(1.8, 3.2, 0, up = FALSE),
  backbone(5.4, 10, 0), ex(5.6, 6.8), ex(7.7, 8.2, event_cols["A3A5"]), ex(8.2, 9.6),
  junction(6.8, 7.7, 0), junction(6.8, 8.2, 0, up = FALSE)))
d_ir <- event_panel("Intron retention (IR)", list(backbone(1, 9, 0), ex(1.5, 3.8), ex(3.8, 6.2, event_cols["IR"]), ex(6.2, 8.5),
  junction(3.8, 6.2, 0, h = 0.4)))
fig1D <- wrap_elements(full = (ggplot() + labs(title = expression(underline("Event-level"))) + schematic_theme) /
                         d_es / d_mx / d_a35 / d_ir + plot_layout(heights = c(0.15, 1, 1, 1, 1)))

# B: isoforms detected by each dataset
fig1B <- readRDS(file.path(fig1_dir, "fig1B_isoforms.rds")) + labs(y = "Isoforms") + fig1_text +
  theme(axis.title.y = element_text(size = 12, angle = 90, margin = margin(r = 8)))

# C: gene expression, long vs short reads; colours as adjusted by hand in the preprint figure
fig1C_cols <- c("#1A151D", "#302840", "#3F3E61", "#485682", "#536E95", "#5C879D", "#6BA0A5", "#7DB9AB", "#A0D0BB", "#CDE7D8")
fig1C <- readRDS(here::here("data", "03_data_analysis", "07_data_correlation", "plots", "fig1C_corr.rds")) +
  scale_fill_gradientn(colours = fig1C_cols, trans = "log10") +
  labs(fill = "Gene\ncount", x = "Long-read\nexpression", y = "Short-read\nexpression") +  # 2-line titles fit the short row
  fig1_text + theme(legend.title = element_text(size = 9), legend.text = element_text(size = 8),
                    # narrow legend so the plot keeps most of the width
                    legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.45, "cm"),
                    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(3, "pt"),
                    plot.margin = margin(16, 2, 2, 22),  # clear of the panel letter
                    aspect.ratio = 0.8)  # a little wider than tall to fit the short row
for (l in fig1C$layers) {  # smaller r label
  if (inherits(l$geom, "GeomText")) l$aes_params$size <- 10 / .pt
}

# E, F: events (and significant events) detected by each dataset, one row per event type
event_stack <- function(file, ytitle) {
  rows <- readRDS(file.path(fig1_dir, file))
  rows <- lapply(rows, function(p) p + labs(y = ytitle) + fig1_text +
                   theme(axis.title.y = element_text(size = 12, angle = 90, face = "plain", colour = "black", margin = margin(r = 8))))
  wrap_plots(rows, ncol = 1) + plot_layout(axis_titles = "collect")
}
fig1E <- wrap_elements(full = event_stack("fig1E_events.rds", "Events"))
fig1F <- wrap_elements(full = event_stack("fig1F_sig_events.rds", "Significant events"))

# One 2 x 3 grid so B sits directly above E and C above F, at the same widths
# (C wrapped so its axis title is not pushed down to line up with B's angled labels)
fig1 <- wrap_elements(full = fig1A) + fig1B + wrap_elements(full = fig1C) +
  fig1D + fig1E + fig1F +
  plot_layout(ncol = 3, widths = c(1, 0.85, 0.85), heights = c(0.33, 1)) +  # top row sized so B matches one E/F row; E/F keep their size
  plot_annotation(tag_levels = "A") &
  tag_theme

ggsave(file.path(out_dir, "14_fig1.png"), plot = fig1, width = 9, height = 6.05, dpi = 300)
ggsave(file.path(out_dir, "14_fig1.pdf"), plot = fig1, width = 9, height = 6.05, device = quartz_pdf)
