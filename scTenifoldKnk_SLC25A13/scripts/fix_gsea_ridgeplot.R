# Fix GSEA ridgeplot - load saved results and regenerate
cat("=== Fixing GSEA Ridgeplot ===\n")

library(data.table)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# Check what GSEA results we saved
gsea_df <- fread("SLC25A13_KO_GSEA.txt", sep = "\t")
cat(sprintf("GSEA results loaded: %d rows\n", nrow(gsea_df)))

# Count significant
sig_gsea <- gsea_df[which(gsea_df$p.adjust < 0.05), ]
cat(sprintf("Significant GSEA gene sets (p.adjust < 0.05): %d\n", nrow(sig_gsea)))

if (nrow(sig_gsea) > 0) {
  cat("Top GSEA terms:\n")
  for (i in 1:min(10, nrow(sig_gsea))) {
    cat(sprintf("  %d. %s (NES=%.3f, p.adjust=%.2e)\n",
                i, sig_gsea$Description[i], sig_gsea$NES[i], sig_gsea$p.adjust[i]))
  }
}

# Reload original KO result to reconstruct GSEA object
result <- readRDS("scTenifoldKnk_SLC25A13_result.rds")
df <- result$diffRegulation
df <- df[df$gene != "SLC25A13", ]

# Recreate ranked gene list
gene_rank <- df$Z
names(gene_rank) <- df$gene
gene_rank <- sort(gene_rank, decreasing = TRUE)

# Re-run GSEA
gsea_all <- gseGO(
  geneList = gene_rank,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  seed = 42,
  verbose = FALSE
)

cat(sprintf("\nGSEA re-run complete: %d gene sets tested\n", nrow(gsea_all@result)))

# Generate ridgeplot - handle NAs safely
if (nrow(gsea_all@result) > 0 && any(!is.na(gsea_all@result$p.adjust) & gsea_all@result$p.adjust < 0.05)) {
  p4 <- ridgeplot(gsea_all, showCategory = min(15, sum(gsea_all@result$p.adjust < 0.05, na.rm = TRUE)))
  ggsave("SLC25A13_KO_GSEA_ridge.pdf", p4, width = 12, height = 10)
  ggsave("SLC25A13_KO_GSEA_ridge.png", p4, width = 12, height = 10, dpi = 150)
  cat("GSEA ridgeplot saved: SLC25A13_KO_GSEA_ridge.pdf/png\n")
} else {
  cat("GSEA ridgeplot: no significant gene sets to plot\n")
}

# Generate GSEA barplot as alternative
if (nrow(gsea_all@result) > 0 && any(!is.na(gsea_all@result$p.adjust) & gsea_all@result$p.adjust < 0.05)) {
  tryCatch({
    p5 <- dotplot(gsea_all, showCategory = min(20, sum(gsea_all@result$p.adjust < 0.05, na.rm = TRUE))) +
      ggtitle("GSEA Biological Process (SLC25A13 KO)")
    ggsave("SLC25A13_KO_GSEA_dotplot.pdf", p5, width = 12, height = 10)
    ggsave("SLC25A13_KO_GSEA_dotplot.png", p5, width = 12, height = 10, dpi = 150)
    cat("GSEA dotplot saved: SLC25A13_KO_GSEA_dotplot.pdf/png\n")
  }, error = function(e) cat("  GSEA dotplot error:", e$message, "\n"))

  # GSEA enrichment plot for top gene set
  tryCatch({
    p6 <- gseaplot2(gsea_all, geneSetID = 1, title = sig_gsea$Description[1])
    ggsave("SLC25A13_KO_GSEA_enrichment_top.pdf", p6, width = 10, height = 7)
    ggsave("SLC25A13_KO_GSEA_enrichment_top.png", p6, width = 10, height = 7, dpi = 150)
    cat("GSEA enrichment plot saved for top gene set\n")
  }, error = function(e) cat("  GSEA enrichment plot error:", e$message, "\n"))
}

cat("\n=== GSEA Ridgeplot Fix Complete ===\n")
