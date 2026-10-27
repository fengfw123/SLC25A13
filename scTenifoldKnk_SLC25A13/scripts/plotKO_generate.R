# KO Subnetwork Visualization
cat("=== KO Subnetwork Visualization ===\n")

library(data.table)
library(ggplot2)
library(ggrepel)

# Load data
df <- fread("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/SLC25A13_KO_all_genes.txt")
cat(sprintf("Total genes: %d\n", nrow(df)))
cat("Column names:", names(df), "\n")

# Filter significant genes
df_sig <- df[df$FC > 2 & df$p.value < 0.05, ]
df_sig <- df_sig[order(df_sig$p.value), ]
top_n <- min(40, nrow(df_sig))
top_genes <- df_sig[1:top_n, ]
cat(sprintf("Top %d significant genes\n", nrow(top_genes)))

# Load GO BP for annotation
go_bp <- fread("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/SLC25A13_KO_GO_BP_301.txt")
top_go <- go_bp[order(go_bp$p.adjust), ][1:min(10, nrow(go_bp)), ]
cat(sprintf("Top GO terms: %d\n", nrow(top_go)))

# Assign GO terms to genes
gene_go_list <- list()
for (i in 1:nrow(top_go)) {
  genes_in_term <- unlist(strsplit(top_go$geneID[i], "/"))
  for (g in genes_in_term) {
    gene_go_list[[length(gene_go_list) + 1]] <- data.table(
      gene = g,
      go_term = top_go$Description[i],
      go_padj = top_go$p.adjust[i]
    )
  }
}
gene_go <- rbindlist(gene_go_list)
gene_go <- unique(gene_go)
# Keep best GO per gene
gene_go <- gene_go[order(gene_go$go_padj), ]
gene_go <- gene_go[!duplicated(gene_go$gene), ]

# ============================================================
# Radial subnetwork plot
# ============================================================
cat("Generating radial subnetwork...\n")

n_plot <- nrow(top_genes)
angle <- seq(0, 2 * pi * (1 - 1/n_plot), length.out = n_plot)
radius <- 8

node_pos <- data.frame(
  id = top_genes$gene,
  x = radius * cos(angle),
  y = radius * sin(angle),
  FC = top_genes$FC
)

# Merge GO annotation
node_pos <- merge(node_pos, gene_go, by.x = "id", by.y = "gene", all.x = TRUE)
node_pos$go_term[is.na(node_pos$go_term)] <- "Other"

# Simplify GO term labels for legend
node_pos$go_simple <- node_pos$go_term
node_pos$go_simple <- gsub("negative regulation of ", "Neg. reg. of ", node_pos$go_simple)
node_pos$go_simple <- gsub("regulation of ", "Reg. of ", node_pos$go_simple)

# Colors
go_cats <- unique(node_pos$go_simple)
set.seed(42)
n_colors <- length(go_cats)
if (n_colors <= 9) {
  pal <- c("#e74c3c","#e67e22","#f1c40f","#2ecc71","#1abc9c",
           "#3498db","#9b59b6","#e84393","#95a5a6")[1:n_colors]
} else {
  pal <- rainbow(n_colors)
}
names(pal) <- go_cats

# Plot
p <- ggplot() +
  # Edges from center
  geom_segment(aes(x = 0, y = 0, xend = node_pos$x, yend = node_pos$y),
               color = "#c0392b", alpha = 0.3, linewidth = 0.4) +
  # SLC25A13 center node
  annotate("point", x = 0, y = 0, size = 16, color = "#302b63", fill = "#302b63", shape = 21) +
  annotate("text", x = 0, y = 0, label = "SLC25A13\nKO", color = "white",
           fontface = "bold", size = 3.5) +
  # Target nodes colored by GO
  geom_point(data = node_pos, aes(x, y, color = go_simple, size = log2(FC)/2),
             alpha = 0.85) +
  # Labels
  geom_text_repel(data = node_pos, aes(x, y, label = id),
                  size = 2.5, max.overlaps = 40, segment.size = 0.2,
                  segment.alpha = 0.5, box.padding = 0.3) +
  scale_color_manual(values = pal, name = "GO BP Category") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  theme_void() +
  coord_fixed() +
  labs(title = "SLC25A13 KO — Gene Regulation Subnetwork",
       subtitle = sprintf("Top %d upregulated genes (FC > 2, p.value < 0.05), colored by GO BP annotation", n_plot),
       caption = "Radial layout: center = SLC25A13 (KO target), edges = significantly upregulated genes") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "#666"),
    legend.position = "bottom",
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_subnetwork.png",
       p, width = 15, height = 13, dpi = 150)
cat("Radial subnetwork PNG saved\n")

ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/pdf/SLC25A13_KO_subnetwork.pdf",
       p, width = 15, height = 13)
cat("Radial subnetwork PDF saved\n")

# ============================================================
# Annotated FC barplot
# ============================================================
cat("Generating annotated FC barplot...\n")

top30 <- top_genes[1:min(30, nrow(top_genes)), ]
top30 <- merge(top30, gene_go, by = "gene", all.x = TRUE)
top30$go_term[is.na(top30$go_term)] <- "Other"
top30$go_simple <- gsub("negative regulation of ", "Neg.reg.", top30$go_term)
top30$go_simple <- gsub("regulation of ", "Reg.of ", top30$go_simple)
top30$go_simple <- gsub("cytoplasmic translation", "Translation", top30$go_simple)

p2 <- ggplot(top30, aes(x = reorder(gene, FC), y = FC, fill = go_simple)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  coord_flip() +
  scale_fill_manual(values = pal, name = "GO BP") +
  labs(title = "Top 30 Upregulated Genes with GO BP Annotation",
       subtitle = "SLC25A13 KO in Melanoma",
       x = "Gene", y = "Fold Change") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_top30_annotated.png",
       p2, width = 12, height = 10, dpi = 150)
cat("Annotated FC barplot saved\n")

ggsave("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/pdf/SLC25A13_KO_top30_annotated.pdf",
       p2, width = 12, height = 10)
cat("Annotated FC barplot PDF saved\n")

# ============================================================
# Try original plotKO via enrichR
# ============================================================
cat("Attempting original scTenifoldKnk::plotKO (via enrichR)...\n")
result <- readRDS("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/scTenifoldKnk_SLC25A13_result.rds")

tryCatch({
  pdf("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/pdf/SLC25A13_KO_enrichr.pdf",
      width = 14, height = 12)
  plotKO(result, gKO = "SLC25A13", q = 0.99, annotate = TRUE,
         nCategories = 30, fdrThreshold = 0.05)
  dev.off()
  cat("plotKO (enrichR) PDF saved\n")

  png("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/images/SLC25A13_KO_enrichr.png",
      width = 1400, height = 1200, res = 150)
  plotKO(result, gKO = "SLC25A13", q = 0.99, annotate = TRUE,
         nCategories = 30, fdrThreshold = 0.05)
  dev.off()
  cat("plotKO (enrichR) PNG saved\n")
}, error = function(e) {
  cat("plotKO (enrichR) failed:", e$message, "\n")
  cat("Manual subnetwork already available instead.\n")
})

cat("\n=== PlotKO Complete ===\n")
