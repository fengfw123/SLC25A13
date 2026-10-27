##############################################################################
# scTenifoldKnk — SLC25A13 KO in Melanoma cells
# Hardware: AMD 5950X (16C/32T) | 128GB RAM | RTX 3070 Ti
# Date: 2026-06-21
##############################################################################

cat("=== scTenifoldKnk: SLC25A13 Knockout in Melanoma ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ============================================================
# 1. Load libraries
# ============================================================
cat("[1/7] Loading libraries...\n")

tryCatch({
  library(scTenifoldKnk)
  cat("  scTenifoldKnk version:", as.character(packageVersion("scTenifoldKnk")), "\n")
}, error = function(e) {
  cat("  Installing scTenifoldKnk from CRAN...\n")
  install.packages("scTenifoldKnk", repos = "https://cloud.r-project.org")
  library(scTenifoldKnk)
})

library(Seurat)
library(data.table)
library(ggplot2)
library(ggrepel)
library(parallel)

# ============================================================
# 2. Load data & subset Melanoma cells
# ============================================================
cat("[2/7] Loading Seurat data...\n")
load("scRNA.umapcell.type(1)+(5).rdata")

obj <- scRNA.norm.umap
rm(scRNA.norm.umap)
gc()

cat(sprintf("  Total cells: %d\n", ncol(obj)))
cat(sprintf("  Total genes: %d\n", nrow(obj)))

# Subset to Melanoma
mel_cells <- WhichCells(obj, ident = "Melanoma")
cat(sprintf("  Melanoma cells: %d\n", length(mel_cells)))

mel_obj <- subset(obj, cells = mel_cells)
rm(obj)
gc()

# ============================================================
# 3. Extract count matrix & filter genes
# ============================================================
cat("[3/7] Extracting count matrix...\n")

counts <- LayerData(mel_obj, assay = "RNA", layer = "counts")
counts <- as.matrix(counts)
cat(sprintf("  Raw matrix: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# Keep top 5000 expressed genes (developer recommended strategy)
gene_means <- rowMeans(counts)
gene_rank <- order(gene_means, decreasing = TRUE)

# Ensure SLC25A13 is included
target_gene <- "SLC25A13"
target_idx <- which(rownames(counts) == target_gene)
cat(sprintf("  '%s' rank by expression: %d / %d\n",
    target_gene, which(gene_rank == target_idx), nrow(counts)))

# Take top 5000, include target if not already in top
top_n <- 5000
top_genes <- rownames(counts)[gene_rank[1:top_n]]
if (!target_gene %in% top_genes) {
  top_genes <- c(target_gene, setdiff(top_genes, target_gene)[1:(top_n-1)])
  cat(sprintf("  Added '%s' back to top %d genes\n", target_gene, top_n))
}

counts_filt <- counts[top_genes, ]
cat(sprintf("  Filtered matrix: %d genes x %d cells\n", nrow(counts_filt), ncol(counts_filt)))
cat(sprintf("  Sparsity: %.1f%% zeros\n",
    sum(counts_filt == 0) / prod(dim(counts_filt)) * 100))

# Free memory
rm(mel_obj)
gc()

# ============================================================
# 4. Run scTenifoldKnk
# ============================================================
cat("[4/7] Running scTenifoldKnk...\n")
cat(sprintf("  Target gene: %s\n", target_gene))
cat(sprintf("  Parameters: nc_nNet=10, nc_nCells=500, nCores=16\n"))
cat(sprintf("  Running... (this may take 2-4 hours)\n\n"))

t_start <- Sys.time()

result <- scTenifoldKnk(
  countMatrix = counts_filt,
  gKO = target_gene,
  qc = FALSE,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3,
  nc_lambda = 0,
  nc_q = 0.9,
  td_K = 3,
  td_maxIter = 1000,
  td_maxError = 1e-05,
  ma_nDim = 2,
  nCores = 16
)

t_end <- Sys.time()
elapsed <- difftime(t_end, t_start, units = "hours")
cat(sprintf("\n  Elapsed time: %.2f hours (%.1f minutes)\n",
    elapsed, as.numeric(elapsed) * 60))

# Save full result
cat("[5/7] Saving results...\n")
saveRDS(result, file = "scTenifoldKnk_SLC25A13_result.rds")

# ============================================================
# 5. Extract differential regulation results
# ============================================================
df <- result$diffRegulation
cat(sprintf("  Total genes evaluated: %d\n", nrow(df)))
cat(sprintf("  Significant (p.adj<0.05): %d\n", sum(df$p.adj < 0.05, na.rm = TRUE)))
cat(sprintf("  Up-regulated (FC>1): %d\n", sum(df$FC > 1 & df$p.adj < 0.05, na.rm = TRUE)))
cat(sprintf("  Down-regulated (FC<1): %d\n", sum(df$FC < 1 & df$p.adj < 0.05, na.rm = TRUE)))

# Remove the KO gene itself from output
df_out <- df[df$gene != target_gene, ]
sig_genes <- df_out[df_out$p.adj < 0.05, ]
sig_genes <- sig_genes[order(sig_genes$FC, decreasing = TRUE), ]

fwrite(df_out, file = "SLC25A13_KO_all_genes.txt", sep = "\t", quote = FALSE)
fwrite(sig_genes, file = "SLC25A13_KO_sig_genes.txt", sep = "\t", quote = FALSE)

# Print top regulated genes
cat("\nTop 20 up-regulated genes (by FC):\n")
top_up <- head(sig_genes[sig_genes$FC > 1, ], 20)
if (nrow(top_up) > 0) {
  print(top_up[, c("gene", "FC", "Z", "p.value", "p.adj")])
}

cat("\nTop 20 down-regulated genes (by FC):\n")
top_down <- head(sig_genes[sig_genes$FC < 1, ], 20)
if (nrow(top_down) > 0) {
  print(top_down[, c("gene", "FC", "Z", "p.value", "p.adj")])
}

# ============================================================
# 6. Visualizations
# ============================================================
cat("[6/7] Generating plots...\n")

# 6a. Volcano plot
df_out$log_p.adj <- -log10(df_out$p.adj)
df_out$significant <- ifelse(df_out$p.adj < 0.05, "Significant", "Not significant")
y_upper <- quantile(df_out$log_p.adj, 0.999, na.rm = TRUE)

# Label top 20 significant genes
sig_label <- df_out[df_out$p.adj < 0.05, ]
sig_label <- sig_label[order(abs(sig_label$Z), decreasing = TRUE), ]
sig_label <- head(sig_label, 20)

p_volcano <- ggplot(df_out, aes(x = Z, y = log_p.adj, color = significant)) +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not significant" = "gray70")) +
  geom_text_repel(data = sig_label, aes(label = gene), size = 3, max.overlaps = 50) +
  labs(title = paste0("SLC25A13 KO - Differential Regulation (Melanoma)"),
       subtitle = paste0(sum(df_out$p.adj < 0.05, na.rm = TRUE), " significant genes"),
       x = "Z-score", y = "-log10(p.adj)") +
  theme_classic() +
  coord_cartesian(ylim = c(0, y_upper)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("SLC25A13_KO_volcano.pdf", p_volcano, width = 8, height = 7)
ggsave("SLC25A13_KO_volcano.png", p_volcano, width = 8, height = 7, dpi = 150)

# 6b. Top 30 FC barplot
top30 <- head(df_out[order(-abs(df_out$FC)), ], 30)
top30$direction <- ifelse(top30$FC > 1, "Up", "Down")
top30$gene <- factor(top30$gene, levels = top30$gene[order(top30$FC)])

p_bar <- ggplot(top30, aes(x = gene, y = FC, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
  labs(title = "Top 30 Genes by Fold Change (SLC25A13 KO)",
       x = "Gene", y = "Fold Change (FC)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

ggsave("SLC25A13_KO_top30_FC.pdf", p_bar, width = 8, height = 7)
ggsave("SLC25A13_KO_top30_FC.png", p_bar, width = 8, height = 7, dpi = 150)

# 6c. Proportion pie
sig_count <- sum(df_out$p.adj < 0.05, na.rm = TRUE)
non_sig_count <- sum(df_out$p.adj >= 0.05, na.rm = TRUE) + sum(is.na(df_out$p.adj))
pie_df <- data.frame(
  category = c("Significant", "Not significant"),
  count = c(sig_count, non_sig_count)
)
pie_df$percentage <- paste0(round(pie_df$count / sum(pie_df$count) * 100, 1), "%")

p_pie <- ggplot(pie_df, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), size = 5) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Significant" = "#E41A1C", "Not significant" = "gray80")) +
  labs(title = "Proportion of Significantly Regulated Genes", fill = "") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right")

ggsave("SLC25A13_KO_proportion.pdf", p_pie, width = 6, height = 5)
ggsave("SLC25A13_KO_proportion.png", p_pie, width = 6, height = 5, dpi = 150)

# ============================================================
# 7. Summary report
# ============================================================
cat("[7/7] Summary report\n\n")
cat("========================================\n")
cat("  scTenifoldKnk — SLC25A13 KO in Melanoma\n")
cat("========================================\n")
cat(sprintf("  Input cells: %d\n", ncol(counts_filt)))
cat(sprintf("  Input genes: %d\n", nrow(counts_filt)))
cat(sprintf("  Execution time: %.2f hours\n", as.numeric(elapsed)))
cat(sprintf("  Significant genes (p.adj<0.05): %d\n", sig_count))
cat(sprintf("  Up-regulated: %d\n", sum(df_out$FC > 1 & df_out$p.adj < 0.05, na.rm = TRUE)))
cat(sprintf("  Down-regulated: %d\n", sum(df_out$FC < 1 & df_out$p.adj < 0.05, na.rm = TRUE)))
cat("========================================\n")
cat("Output files:\n")
cat("  scTenifoldKnk_SLC25A13_result.rds  — full R object\n")
cat("  SLC25A13_KO_all_genes.txt           — all genes with stats\n")
cat("  SLC25A13_KO_sig_genes.txt           — significant genes only\n")
cat("  SLC25A13_KO_volcano.pdf/png         — volcano plot\n")
cat("  SLC25A13_KO_top30_FC.pdf/png        — FC barplot\n")
cat("  SLC25A13_KO_proportion.pdf/png      — proportion pie\n")
cat("========================================\n")

cat("\nEnd time:", format(Sys.time()), "\n")
cat("=== Analysis Complete ===\n")
