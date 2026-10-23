suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("code/functions/common.R")
root <- "results/gse162170/within_class"
x <- fread(file.path(root, "pooled_within_class_overlap.tsv"))
loci <- fread(file.path(root, "focal_locus_overlap.tsv"))
classes <- c("GluN", "IN", "IPC_RG", "OPC_Glia")
labels <- c("GluN\n5,564 cells", "IN\n2,040 cells", "IPC/RG\n987 cells", "OPC/Glia\n359 cells")
x[, cell_class := factor(cell_class, levels = classes)]
loci[, cell_class := factor(cell_class, levels = classes)]
x[, endpoint := factor(endpoint, levels = c("TF-gene", "factorized-chain"), labels = c("Direct TF–gene", "Exact-region chains"))]
a <- ggplot(x, aes(cell_class, signed_jaccard, colour = endpoint, group = endpoint)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  scale_colour_manual(values = c("#2177B8", MULTICSN_COLOUR)) +
  scale_x_discrete(labels = labels) + labs(x = NULL, y = "Signed Jaccard with pooled fit", colour = NULL, tag = "a") +
  theme_genome(7) + theme(legend.position = "top", axis.text.x = element_text(size = 6.3))
loci[, fraction_label := paste0(shared_signed, "/", pooled_chains)]
loci[, label_y := pooled_signed_recall + fifelse(target == "KCNC2", 0.014, 0.034)]
b <- ggplot(loci, aes(cell_class, pooled_signed_recall, colour = target, group = target)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  geom_text(aes(y = label_y, label = fraction_label),
            size = 6 / ggplot2::.pt, show.legend = FALSE) +
  scale_colour_manual(values = c("#8064A2", "#2177B8")) +
  scale_x_discrete(labels = labels) + scale_y_continuous(expand = expansion(mult = c(0.08, 0.2))) +
  labs(x = NULL, y = "Same-sign recovery of pooled locus chains", colour = NULL, tag = "b") +
  theme_genome(7) + theme(legend.position = "top", legend.text = element_text(face = "italic"), axis.text.x = element_text(size = 6.3))
figure <- (a | b) & theme(plot.tag = element_text(face = "bold", size = 9))
save_pdf(figure, "figures/figs12.pdf", 7.2, 2.8)
