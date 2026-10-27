##############################################################################
# Re-enrichment with relaxed threshold: FC>2 & p.value<0.05 (301 genes)
# Show top 15 pathways regardless of significance
##############################################################################

cat("=== Re-enrichment with 301 genes (FC>2, p.value<0.05) ===\n\n")

library(data.table)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

# Load data
df <- fread("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/SLC25A13_KO_all_genes.txt")

# Filter: FC > 2 & p.value < 0.05
sig_genes <- df[FC > 2 & p.value < 0.05]
cat(sprintf("筛选基因: %d\n", nrow(sig_genes)))

# Confirm no down genes meet threshold
n_down <- nrow(df[FC < 0.5 & p.value < 0.05])
cat(sprintf("FC<0.5 & p.value<0.05: %d (无下调基因)\n\n", n_down))

# Save new gene list
fwrite(sig_genes, "D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/SLC25A13_KO_sig_genes_301.txt",
  sep = "\t", quote = FALSE)

# ENTREZ ID mapping
gene_list <- sig_genes$gene
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(sprintf("ENTREZ ID 映射: %d/%d\n", nrow(entrez_ids), length(gene_list)))

# ============================================================
# GO Enrichment
# ============================================================
run_enrich <- function(ont_name) {
  result <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont_name,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500,
    readable = TRUE
  )
  return(result)
}

go_bp <- run_enrich("BP")
go_mf <- run_enrich("MF")
go_cc <- run_enrich("CC")

# KEGG
kegg <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500
)

# ============================================================
# Save all results (all terms)
# ============================================================
save_table <- function(obj, filename) {
  if (!is.null(obj)) {
    fwrite(as.data.frame(obj@result),
      paste0("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/", filename),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

save_table(go_bp, "SLC25A13_KO_GO_BP_301.txt")
save_table(go_mf, "SLC25A13_KO_GO_MF_301.txt")
save_table(go_cc, "SLC25A13_KO_GO_CC_301.txt")
save_table(kegg, "SLC25A13_KO_KEGG_301.txt")

# ============================================================
# Print summary (show top 15, or all if <15)
# ============================================================
print_top15 <- function(obj, label) {
  if (is.null(obj)) { cat(sprintf("\n%s: NULL\n", label)); return() }
  n_sig <- sum(obj@result$p.adjust < 0.05, na.rm = TRUE)
  cat(sprintf("\n======== %s ========\n", label))
  cat(sprintf("总项数: %d | p.adjust<0.05: %d\n", nrow(obj@result), n_sig))

  n_show <- min(15, nrow(obj@result))
  top <- obj@result[order(obj@result$p.adjust), ][1:n_show, ]
  cat(sprintf("展示: Top %d (按 p.adjust 排序)\n\n", n_show))
  for (i in 1:nrow(top)) {
    sig_mark <- ifelse(top$p.adjust[i] < 0.05, " *", "")
    cat(sprintf("  %2d. %s%s\n", i, top$Description[i], sig_mark))
    cat(sprintf("      p.adjust=%.2e, 基因数=%d, 富集倍数=%.1f\n",
        top$p.adjust[i], top$Count[i], top$FoldEnrichment[i]))
  }
}

print_top15(go_bp, "GO Biological Process (301 genes)")
print_top15(go_mf, "GO Molecular Function (301 genes)")
print_top15(go_cc, "GO Cellular Component (301 genes)")
print_top15(kegg, "KEGG Pathway (301 genes)")

# ============================================================
# Generate plots - show top 15 (or all if <15)
# ============================================================
n_bp <- nrow(go_bp@result)
n_cc <- nrow(go_cc@result)
n_kegg <- nrow(kegg@result)

cat("\n\n======== 生成图表 ========\n")

# GO BP
if (n_bp > 0) {
  p1 <- dotplot(go_bp, showCategory = min(15, n_bp)) +
    ggtitle("GO Biological Process (SLC25A13 KO, 301 genes)")
  ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_GO_BP_dotplot_301.png",
    p1, width = 11, height = max(6, min(15, n_bp) * 0.5), dpi = 150)
  cat(sprintf("  GO BP dotplot: %d terms\n", min(15, n_bp)))

  # CNet plot (if >= 5 terms)
  if (n_bp >= 5) {
    tryCatch({
      p1c <- cnetplot(go_bp, showCategory = min(10, n_bp),
                      cex_label_category = 0.8, cex_label_gene = 0.5)
      ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_GO_BP_network_301.png",
        p1c, width = 16, height = 14, dpi = 150)
      cat("  GO BP network plot saved\n")
    }, error = function(e) cat("  cnetplot error:", e$message, "\n"))
  }
}

# GO CC
if (n_cc > 0) {
  p2 <- dotplot(go_cc, showCategory = min(15, n_cc)) +
    ggtitle("GO Cellular Component (SLC25A13 KO, 301 genes)")
  ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_GO_CC_dotplot_301.png",
    p2, width = 11, height = max(6, min(15, n_cc) * 0.5), dpi = 150)
  cat(sprintf("  GO CC dotplot: %d terms\n", min(15, n_cc)))
}

# KEGG
if (n_kegg > 0) {
  p3 <- dotplot(kegg, showCategory = min(15, n_kegg)) +
    ggtitle("KEGG Pathway (SLC25A13 KO, 301 genes)")
  ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_KEGG_dotplot_301.png",
    p3, width = 11, height = max(6, min(15, n_kegg) * 0.5), dpi = 150)
  cat(sprintf("  KEGG dotplot: %d terms\n", min(15, n_kegg)))
}

# ============================================================
# GSEA with full ranking
# ============================================================
cat("\n--- GSEA (全部 4999 基因) ---\n")
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
  pvalueCutoff = 1,
  seed = 42,
  verbose = FALSE
)

save_table(gsea_all, "SLC25A13_KO_GSEA_301.txt")

n_gsea <- nrow(gsea_all@result)
n_sig_gsea <- sum(gsea_all@result$p.adjust < 0.05, na.rm = TRUE)
cat(sprintf("GSEA 总项: %d | p.adjust<0.05: %d\n", n_gsea, n_sig_gsea))

# Show top 15
n_show_gsea <- min(15, n_gsea)
top_gsea <- gsea_all@result[order(gsea_all@result$p.adjust), ][1:n_show_gsea, ]
cat(sprintf("Top %d GSEA terms:\n", n_show_gsea))
for (i in 1:nrow(top_gsea)) {
  sig_mark <- ifelse(top_gsea$p.adjust[i] < 0.05, " *", "")
  cat(sprintf("  %2d. %s%s (NES=%.2f, p.adj=%.2e)\n",
      i, top_gsea$Description[i], sig_mark, top_gsea$NES[i], top_gsea$p.adjust[i]))
}

# GSEA plots
if (n_gsea > 0 && any(!is.na(gsea_all@result$p.adjust))) {
  tryCatch({
    p4 <- ridgeplot(gsea_all, showCategory = min(15, n_gsea)) +
      ggtitle("GSEA Biological Process (SLC25A13 KO)")
    ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_GSEA_ridge_301.png",
      p4, width = 14, height = max(7, min(15, n_gsea) * 0.5), dpi = 150)
    cat("  GSEA ridgeplot saved\n")
  }, error = function(e) cat("  ridgeplot error:", e$message, "\n"))
}

cat("\n=== Re-enrichment Complete ===\n")
