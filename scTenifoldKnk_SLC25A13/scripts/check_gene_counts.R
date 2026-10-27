# 检查不同阈值的基因数量
df <- read.table("D:/03_Agent/Paper_xjn/SLC25A13_KO_Report/data/SLC25A13_KO_all_genes.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("总基因数: %d\n", nrow(df)))
cat(sprintf("SLC25A13 已排除? %s\n", ifelse(!"SLC25A13" %in% df$gene, "是", "否")))

# 条件1: 原阈值 p.adj < 0.05
n1 <- sum(df$p.adj < 0.05, na.rm = TRUE)
cat(sprintf("\n--- 原阈值 p.adj < 0.05 ---\n"))
cat(sprintf("  显著基因: %d\n", n1))
cat(sprintf("  其中 FC>1 (上调): %d\n", sum(df$p.adj < 0.05 & df$FC > 1, na.rm = TRUE)))
cat(sprintf("  其中 FC<1 (下调): %d\n", sum(df$p.adj < 0.05 & df$FC < 1, na.rm = TRUE)))

# 条件2: FC > 2 且 p.adj < 0.05
n2_up <- sum(df$FC > 2 & df$p.adj < 0.05, na.rm = TRUE)
n2_down <- sum(df$FC > 2 & df$p.value < 0.05 & df$p.adj >= 0.05, na.rm = TRUE)
cat(sprintf("\n--- FC > 2 ---\n"))
cat(sprintf("  + p.adj < 0.05: %d 个\n", n2_up))
# Check FC for down-regulated
n2_down_adj <- sum(df$FC < 0.5 & df$p.adj < 0.05, na.rm = TRUE)
n2_down_p <- sum(df$FC < 0.5 & df$p.value < 0.05 & df$p.adj >= 0.05, na.rm = TRUE)
cat(sprintf("  FC < 0.5 + p.adj < 0.05: %d 个\n", n2_down_adj))
cat(sprintf("  FC < 0.5 + p.value < 0.05 (不含上述): %d 个\n", n2_down_p))

# 条件3: 宽松阈值 FC>2 且 (p.adj<0.05 | p.value<0.05)
n3_up <- sum(df$FC > 2 & (df$p.adj < 0.05 | df$p.value < 0.05), na.rm = TRUE)
n3_down <- sum(df$FC < 0.5 & (df$p.adj < 0.05 | df$p.value < 0.05), na.rm = TRUE)
cat(sprintf("\n--- FC>2 | FC<0.5 且 (p.adj<0.05 | p.value<0.05) ---\n"))
cat(sprintf("  上调 (FC>2): %d\n", n3_up))
cat(sprintf("  下调 (FC<0.5): %d\n", n3_down))
cat(sprintf("  总计: %d\n", n3_up + n3_down))

# 看看是否有 FC<1 的基因
n_under1 <- sum(df$FC < 1, na.rm = TRUE)
n_under05 <- sum(df$FC < 0.5, na.rm = TRUE)
cat(sprintf("\n--- FC 分布 ---\n"))
cat(sprintf("  FC < 1: %d 个\n", n_under1))
cat(sprintf("  FC < 0.5: %d 个\n", n_under05))
cat(sprintf("  FC < 0.8: %d 个\n", sum(df$FC < 0.8, na.rm = TRUE)))
cat(sprintf("  最小 FC: %.4f\n", min(df$FC, na.rm = TRUE)))
cat(sprintf("  最大 FC: %.4f\n", max(df$FC, na.rm = TRUE)))
cat(sprintf("  FC 中位数: %.4f\n", median(df$FC, na.rm = TRUE)))
cat(sprintf("  FC < 1 的基因 p.adj < 0.05: %d\n", sum(df$FC < 1 & df$p.adj < 0.05, na.rm = TRUE)))

# 列出所有 FC<1 的基因
down_genes <- df[df$FC < 1, c("gene", "FC", "p.value", "p.adj")]
if (nrow(down_genes) > 0) {
  cat(sprintf("\n--- 所有 FC<1 的基因 (n=%d) ---\n", nrow(down_genes)))
  for (i in 1:nrow(down_genes)) {
    cat(sprintf("  %s  FC=%.4f  p.value=%.2e  p.adj=%.2e\n",
        down_genes$gene[i], down_genes$FC[i], down_genes$p.value[i], down_genes$p.adj[i]))
  }
}

cat(sprintf("\n\n推荐阈值 (FC>2 | FC<0.5) & (p.adj<0.05 | p.value<0.05):\n"))
cat(sprintf("  上调: %d, 下调: %d, 合计: %d\n", n3_up, n3_down, n3_up + n3_down))
