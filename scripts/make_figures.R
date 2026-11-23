# Three figures for the MAGE comparison of T2T-CHM13 and GRCh38.
#
#   Rscript scripts/make_figures.R <bundle> <paired-test dir> <gene-sets table> <output dir>
#
# Called by the pipeline process FIGURES. Inputs:
#   <bundle>/{star,salmon}/egenes.tsv           eGene status per reference (compare_eqtl.py, via export_bundle.py)
#   <paired-test dir>/effects-star.tsv          same-variant effects on both references (paired_test.py)
#   <gene-sets table>                           changed genes and their tier (gene_sets.py)
# Writes <output dir>/fig0N-*.{png,pdf}:
#   fig01-egenes           eGenes (FDR 5%) per reference and quantification method
#   fig02-effects          eQTL effect of the same variant on GRCh38 (x) and T2T (y); changed genes highlighted
#   fig03-changed-genes    changed genes (main set) by region class
# Style: journal figures (<= 183 mm wide, 7 pt text, no titles), Okabe-Ito colours.
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(scales) })
args <- commandArgs(TRUE)
if (length(args) != 4) stop("usage: make_figures.R <bundle> <paired-test dir> <gene-sets table> <output dir>")
bundle <- args[1]; paired <- args[2]; gene_sets <- args[3]; fd <- args[4]
dir.create(fd, showWarnings = FALSE, recursive = TRUE)
rd <- function(path) fread(path, sep = "\t", colClasses = list(character = "gene_ncbi"))

ref_lab <- c("chm13v2.0" = "T2T-CHM13", "grch38" = "GRCh38")
ref_col <- c("T2T-CHM13" = "#0072B2", "GRCh38" = "#D55E00")
reg_lab <- c(other = "Other", sd = "Segmental\nduplication", t2t_new_sequence = "T2T-only\nsequence")

# TeX Gyre Heros (Helvetica metrics; scripts/install_fonts.sh) if installed, else the default sans font
FONT <- if (any(grepl("TeX Gyre Heros", tryCatch(system("fc-list : family", intern = TRUE), error = function(e) ""))))
  "TeX Gyre Heros" else "sans"
pdf(NULL)  # no stray Rplots.pdf
nat <- theme_classic(base_size = 7, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        axis.ticks.length = unit(1.2, "pt"), axis.text = element_text(size = 6.5, colour = "black"),
        axis.title = element_text(size = 7), legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        legend.key.size = unit(7, "pt"), legend.background = element_blank(), plot.margin = margin(2, 4, 2, 2))
save_nat <- function(p, name, w_mm, h_mm) {
  ggsave(file.path(fd, paste0(name, ".png")), p, width = w_mm, height = h_mm, units = "mm", dpi = 300, bg = "white",
         device = png, type = "cairo")
  ggsave(file.path(fd, paste0(name, ".pdf")), p, width = w_mm, height = h_mm, units = "mm", device = cairo_pdf)
}

# Fig 1: eGenes per reference ---------------------------------------------------------------------------------------
eg <- rbindlist(lapply(c(star = "STAR", salmon = "salmon"), function(m) {
  g <- rd(file.path(bundle, tolower(m), "egenes.tsv"))
  data.table(method = m, ref = ref_lab, egenes = c(sum(g$egene_chm13v2.0, na.rm = TRUE), sum(g$egene_grch38, na.rm = TRUE)))
}))
eg[, method := factor(method, levels = c("STAR", "salmon"))]
eg[, ref := factor(ref, levels = ref_lab)]
p1 <- ggplot(eg, aes(method, egenes, fill = ref)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_text(aes(label = comma(egenes)), position = position_dodge(0.75), vjust = -0.4, size = 6 / .pt, family = FONT) +
  scale_fill_manual(values = ref_col, name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "eGenes (FDR 5%)") + nat + theme(legend.position = "bottom")
save_nat(p1, "fig01-egenes", 70, 60)

# Fig 2: same variant, two references -------------------------------------------------------------------------------
ef <- rd(file.path(paired, "effects-star.tsv"))[artefact_flag == FALSE]
gs <- rd(gene_sets)
ef[, grp := fifelse(gene_ncbi %in% gs[tier %in% c("core", "main"), gene_ncbi], "Changed", "Not changed")]
n2 <- ef[, .N, by = grp]
lab2 <- setNames(sprintf("%s (%s genes)", n2$grp, comma(n2$N)), n2$grp)
lim <- max(abs(c(ef$slope_a, ef$slope_b)), na.rm = TRUE)
p2 <- ggplot(ef[order(grp == "Changed")], aes(slope_b, slope_a, colour = grp, size = grp)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") + geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
  geom_point(stroke = 0) +
  scale_size_manual(values = c("Changed" = 1.1, "Not changed" = 0.35), guide = "none") +
  scale_colour_manual(values = c("Changed" = "#0072B2", "Not changed" = "grey80"), labels = lab2, name = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.6))) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "Effect on GRCh38 (SD per allele)", y = "Effect on T2T-CHM13 (SD per allele)") +
  nat + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
save_nat(p2, "fig02-effects", 85, 85)

# Fig 3: changed genes by region class ------------------------------------------------------------------------------
ch <- gs[tier %in% c("core", "main")]
ch[, region := factor(fifelse(region_class %in% names(reg_lab), reg_lab[region_class], region_class),
                      levels = unique(c(reg_lab, region_class)))]
c3 <- ch[, .N, by = region]
p3 <- ggplot(c3, aes(region, N)) + geom_col(width = 0.6, fill = "#0072B2") +
  geom_text(aes(label = N), vjust = -0.4, size = 6.5 / .pt, family = FONT) +
  scale_x_discrete(drop = TRUE) + scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Changed genes, T2T vs GRCh38") + nat
save_nat(p3, "fig03-changed-genes", 70, 60)

message(sprintf("eGenes %s; %d genes compared, %d changed; figures in %s",
                paste(sprintf("%s %s %s", eg$method, eg$ref, comma(eg$egenes)), collapse = ", "),
                nrow(ef), nrow(ch), fd))
