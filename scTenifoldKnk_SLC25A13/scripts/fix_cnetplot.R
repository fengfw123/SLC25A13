# Fix cnetplot with updated enrichplot API
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

# Load data
df <- fread("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/SLC25A13_KO_all_genes.txt")
sig_genes <- df[FC > 2 & p.value < 0.05]
gene_list <- sig_genes$gene
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Re-run GO BP
go_bp <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 500,
  readable = TRUE
)

# cnetplot with updated parameters
n_bp <- nrow(go_bp@result)
cat(sprintf("GO BP terms: %d\n", n_bp))

if (n_bp >= 5) {
  tryCatch({
    p <- cnetplot(go_bp, showCategory = min(10, n_bp),
                  size_category = 0.8, size_item = 0.5)
    ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_GO_BP_network_301.png",
      p, width = 16, height = 14, dpi = 150)
    cat("cnetplot saved successfully\n")
  }, error = function(e) {
    cat("cnetplot error:", e$message, "\n")
    # fallback: try with defaults
    tryCatch({
      p <- cnetplot(go_bp, showCategory = min(10, n_bp))
      ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_GO_BP_network_301.png",
        p, width = 16, height = 14, dpi = 150)
      cat("cnetplot (fallback) saved\n")
    }, error = function(e2) cat("cnetplot fallback also failed:", e2$message, "\n"))
  })
}
