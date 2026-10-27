##############################################################################
# Deep Analysis of SLC25A13 KO Results
# 1. GO/KEGG enrichment (clusterProfiler)
# 2. KO subnetwork visualization
# 3. Integration with melanoma biology
# 4. GSEA analysis of full gene ranking
##############################################################################

cat("=== Post-KO Deep Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ============================================================
# 1. Load libraries
# ============================================================
cat("[1/7] Loading libraries...\n")

library(Seurat)
library(data.table)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)

# Check clusterProfiler
if (!require(clusterProfiler)) {
  if (!require(BiocManager)) install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
}
if (!require(enrichplot)) {
  BiocManager::install("enrichplot")
  library(enrichplot)
}
library(scTenifoldKnk)

# ============================================================
# 2. Load KO results
# ============================================================
cat("[2/7] Loading KO results...\n")

result <- readRDS("scTenifoldKnk_SLC25A13_result.rds")

df <- result$diffRegulation
df <- df[df$gene != "SLC25A13", ]

sig_genes <- df[df$p.adj < 0.05, ]
sig_up <- sig_genes[sig_genes$FC > 1, ]
sig_down <- sig_genes[sig_genes$FC < 1, ]

cat(sprintf("  Significant genes: %d\n", nrow(sig_genes)))
cat(sprintf("  Up-regulated: %d\n", nrow(sig_up)))
cat(sprintf("  Down-regulated: %d\n", nrow(sig_down)))

# ============================================================
# 3. GO Enrichment Analysis (BP, CC, MF, KEGG)
# ============================================================
cat("[3/7] Running GO enrichment...\n")

# Convert gene symbols to ENTREZ ID
gene_list <- sig_genes$gene
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(sprintf("  Mapped %d/%d genes to ENTREZ IDs\n", nrow(entrez_ids), length(gene_list)))

# GO Biological Process
go_bp <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)

# GO Molecular Function
go_mf <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)

# GO Cellular Component
go_cc <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)

# KEGG Pathway
kegg <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 5,
  maxGSSize = 500
)

# Save enrichment results
if (!is.null(go_bp) && nrow(go_bp@result[go_bp@result$p.adjust < 0.05, ]) > 0) {
  fwrite(as.data.frame(go_bp@result), "SLC25A13_KO_GO_BP.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  GO BP: %d significant terms\n", sum(go_bp@result$p.adjust < 0.05)))
} else {
  cat("  GO BP: no significant terms\n")
}

if (!is.null(go_mf) && nrow(go_mf@result[go_mf@result$p.adjust < 0.05, ]) > 0) {
  fwrite(as.data.frame(go_mf@result), "SLC25A13_KO_GO_MF.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  GO MF: %d significant terms\n", sum(go_mf@result$p.adjust < 0.05)))
}

if (!is.null(go_cc) && nrow(go_cc@result[go_cc@result$p.adjust < 0.05, ]) > 0) {
  fwrite(as.data.frame(go_cc@result), "SLC25A13_KO_GO_CC.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  GO CC: %d significant terms\n", sum(go_cc@result$p.adjust < 0.05)))
}

if (!is.null(kegg) && nrow(kegg@result[kegg@result$p.adjust < 0.05, ]) > 0) {
  fwrite(as.data.frame(kegg@result), "SLC25A13_KO_KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  KEGG: %d significant pathways\n", sum(kegg@result$p.adjust < 0.05)))
}

# ============================================================
# 4. GSEA (using full gene ranking by Z-score)
# ============================================================
cat("[4/7] Running GSEA...\n")

# Create ranked gene list (signed by FC direction * Z magnitude)
gene_rank <- df$Z
names(gene_rank) <- df$gene
gene_rank <- sort(gene_rank, decreasing = TRUE)

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

if (!is.null(gsea_all) && nrow(gsea_all@result[gsea_all@result$p.adjust < 0.05, ]) > 0) {
  fwrite(as.data.frame(gsea_all@result), "SLC25A13_KO_GSEA.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  GSEA: %d significant gene sets\n", sum(gsea_all@result$p.adjust < 0.05)))
} else {
  cat("  GSEA: no significant gene sets\n")
}

# ============================================================
# 5. KO Subnetwork Visualization
# ============================================================
cat("[5/7] Generating KO subnetwork plot...\n")

# Plot KO-centered subnetwork using scTenifoldKnk built-in function
tryCatch({
  pdf("SLC25A13_KO_subnetwork.pdf", width = 12, height = 10)
  plotKO(result, gKO = "SLC25A13", top = 30)
  dev.off()
  cat("  KO subnetwork saved: SLC25A13_KO_subnetwork.pdf\n")
}, error = function(e) {
  cat("  plotKO not available (may need different scTenifoldKnk version)\n")
})

# ============================================================
# 6. Enrichment Visualization Plots
# ============================================================
cat("[6/7] Generating enrichment plots...\n")

# GO BP dotplot
if (exists("go_bp") && !is.null(go_bp) && sum(go_bp@result$p.adjust < 0.05) > 0) {
  p1 <- dotplot(go_bp, showCategory = 20, title = "GO Biological Process (SLC25A13 KO)")
  ggsave("SLC25A13_KO_GO_BP_dotplot.pdf", p1, width = 10, height = 8)
  ggsave("SLC25A13_KO_GO_BP_dotplot.png", p1, width = 10, height = 8, dpi = 150)
}

# GO BP network plot
if (exists("go_bp") && !is.null(go_bp) && sum(go_bp@result$p.adjust < 0.05) > 5) {
  tryCatch({
    p2 <- cnetplot(go_bp, showCategory = 10,
                   cex_label_category = 0.8, cex_label_gene = 0.6)
    ggsave("SLC25A13_KO_GO_BP_network.pdf", p2, width = 14, height = 12)
    ggsave("SLC25A13_KO_GO_BP_network.png", p2, width = 14, height = 12, dpi = 150)
  }, error = function(e) cat("  cnetplot skipped\n"))
}

# KEGG dotplot
if (exists("kegg") && !is.null(kegg) && sum(kegg@result$p.adjust < 0.05) > 0) {
  p3 <- dotplot(kegg, showCategory = 15, title = "KEGG Pathways (SLC25A13 KO)")
  ggsave("SLC25A13_KO_KEGG_dotplot.pdf", p3, width = 10, height = 7)
  ggsave("SLC25A13_KO_KEGG_dotplot.png", p3, width = 10, height = 7, dpi = 150)
}

# GSEA ridgeplot
if (exists("gsea_all") && !is.null(gsea_all) && sum(gsea_all@result$p.adjust < 0.05) > 0) {
  p4 <- ridgeplot(gsea_all, showCategory = 15)
  ggsave("SLC25A13_KO_GSEA_ridge.pdf", p4, width = 12, height = 10)
  ggsave("SLC25A13_KO_GSEA_ridge.png", p4, width = 12, height = 10, dpi = 150)
}

# ============================================================
# 7. Summary Report
# ============================================================
cat("[7/7] Summary report\n\n")

cat("========================================\n")
cat("  SLC25A13 KO — Deep Analysis Summary\n")
cat("========================================\n")

# Print GO BP top terms
if (exists("go_bp") && !is.null(go_bp) && sum(go_bp@result$p.adjust < 0.05) > 0) {
  top_bp <- head(go_bp@result[go_bp@result$p.adjust < 0.05, ], 10)
  cat("\nTop GO Biological Process terms:\n")
  for (i in 1:nrow(top_bp)) {
    cat(sprintf("  %d. %s (p.adjust=%.2e, %d genes)\n",
        i, top_bp$Description[i], top_bp$p.adjust[i], top_bp$Count[i]))
  }
}

if (exists("kegg") && !is.null(kegg) && sum(kegg@result$p.adjust < 0.05) > 0) {
  top_kegg <- head(kegg@result[kegg@result$p.adjust < 0.05, ], 10)
  cat("\nTop KEGG Pathways:\n")
  for (i in 1:nrow(top_kegg)) {
    cat(sprintf("  %d. %s (p.adjust=%.2e, %d genes)\n",
        i, top_kegg$Description[i], top_kegg$p.adjust[i], top_kegg$Count[i]))
  }
}

cat("\n========================================\n")
cat("Output files:\n")
cat("  SLC25A13_KO_GO_BP.txt          — GO Biological Process enrichment\n")
cat("  SLC25A13_KO_GO_MF.txt          — GO Molecular Function enrichment\n")
cat("  SLC25A13_KO_GO_CC.txt          — GO Cellular Component enrichment\n")
cat("  SLC25A13_KO_KEGG.txt           — KEGG pathway enrichment\n")
cat("  SLC25A13_KO_GSEA.txt           — GSEA results (full gene ranking)\n")
cat("  SLC25A13_KO_GO_BP_dotplot.*    — GO dotplot\n")
cat("  SLC25A13_KO_GO_BP_network.*    — GO gene-concept network\n")
cat("  SLC25A13_KO_KEGG_dotplot.*     — KEGG dotplot\n")
cat("  SLC25A13_KO_GSEA_ridge.*       — GSEA ridgeplot\n")
cat("  SLC25A13_KO_subnetwork.pdf     — KO-centered subnetwork\n")
cat("========================================\n")

cat("\nEnd time:", format(Sys.time()), "\n")
cat("=== Deep Analysis Complete ===\n")
