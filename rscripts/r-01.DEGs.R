rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "01_DEGs"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)

pacman::p_load(limma, magrittr, stringr, tidyverse, readr, 
               psych, ggplot2, ggthemes, RColorBrewer, Ipaper, scales ,ggrepel,
               ComplexHeatmap)
############################       DEGs      ###################################
# 01.原始数据 ------------------------------------------------------------------
gse <- "GSE15605"
diseases <- "SKCM"
dat <- read.csv(paste0('../00_rawdata/',gse,'/01.dat(',gse,').csv'),row.names = 1)
group <- read.csv(paste0('../00_rawdata/',gse,'/02.group(',gse,').csv'))
dat <- dat[,group$sample]
table(group$group)
# Control    SKCM 
# 16      28 

group$group <- factor(group$group, levels = c("Control", "SKCM"))
design.mat <- cbind(Control = ifelse(group$group == "Control", 1, 0), 
                    SKCM = ifelse(group$group == "Control", 0, 1))

#  差异分析
contrast.mat <- makeContrasts(contrasts="SKCM-Control", levels=design.mat)
fit <- lmFit(dat, design.mat)
fit <- contrasts.fit(fit, contrast.mat)
fit <- eBayes(fit)
fit <- topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")

DEG <- na.omit(fit)


### 根据阈值筛选基因
### |log2FC| > 1.5，P < 0.05
logFC_cutoff <- 1
colnames(DEG)
DEG$change<- as.factor(
  ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
#  19492           
nrow(sig_diff)
paste0("最终获得 ",nrow(sig_diff)," 个差异表达基因")
### "最终获得 2826 个差异表达基因"
table(sig_diff$change)
### DOWN  NOT   UP 
### 1689    0 1137    

write.csv(DEG,file = paste0("01.DEGs_all(",gse,").csv"))
write.csv(sig_diff,file = paste0("02.DEGs_sig(",gse,"2).csv"))


#  02.火山图 -------------------------------------------------------------------
dat_rep <- DEG[rownames(DEG)%in%
                 rownames(rbind(head(sig_diff[order(sig_diff$logFC,
                                                    decreasing = T),],10),
                                head(sig_diff[order(sig_diff$logFC,
                                                    decreasing = F),],10))),]
DEG$Symbols <- rownames(DEG)
up_DEG <- DEG[which(DEG$change == "UP"),]
up_DEG <- up_DEG[order(-up_DEG$logFC, up_DEG$P.Value),]
down_DEG <- DEG[which(DEG$change == "DOWN"),]
down_DEG <- down_DEG[order(down_DEG$logFC, down_DEG$P.Value),]
# data_repel <- rbind(up_DEG[1:10,], down_DEG[1:10,])

DEG$change <- factor(DEG$change, levels =  c("UP", "NOT", "DOWN"))
max_lfc <- max(DEG$logFC)
min_lfc <- min(DEG$logFC)
# max_adj.P.Val <- round(max(-log10(DEG$pvalue)))

if(T){
  volcano_plot <- ggplot(data = DEG, 
                         aes(x = logFC,
                             y = -log10(P.Value), 
                             color =change)) +
    scale_color_manual(values = c("#ed556a", "darkgray","#104680")) +
    geom_point(size = 2, alpha = 0.6, na.rm=T) +
    theme_bw(base_size = 12) +
    geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff),
               lty = 4,col = "darkgray",lwd = 0.6)+
    geom_hline(yintercept = -log10(0.05),
               lty = 4,col = "darkgray",lwd = 0.6)+
    theme(legend.position = "right",
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face="bold",color="black",size=13),
          plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 18),
          axis.text.x = element_text(face = "bold",color = "black",size = 15),
          axis.text.y = element_text(face = "bold",color = "black",size = 15),
          axis.title.x = element_text(face = "bold",color = "black",size = 15),
          axis.title.y = element_text(face = "bold",color = "black",size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic", colour = "black")) +
    geom_label_repel(
      data = dat_rep,
      aes(label = rownames(dat_rep)),
      max.overlaps = 20,
      size = 4,
      box.padding = unit(0.5, "lines"),
      min.segment.length = 0,
      point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
    labs(x = "log2(Fold Change)",y = "-log10 (P.Value)",  ### P.Value ， adj.P.Val
         title = paste0(diseases," vs Control"),
         subtitle = paste(sprintf('P.Value: %.2f;', 0.05), 
                          sprintf('log2 (FoldChange): %.2f;', logFC_cutoff),
                          sprintf('Up: %1.0f; Down: %1.0f;', dim(up_DEG)[1], dim(down_DEG)[1]),
                          sprintf('Total: %1.0f', dim(sig_diff)[1]))) 
  volcano_plot
  ggsave(paste0('03.volcano(',gse,').pdf'), volcano_plot,width = 8, height = 6)
  ggsave(paste0('03.volcano(',gse,').png'), volcano_plot,width = 8, height = 6, dpi = 600, units = "in")
}


# 03.热图  ---------------------------------------------------------------------
rt <- dat
rt <- rt[rownames(dat_rep),]
# rt <- rt[rownames(sig_diff),]
group <- group[order(group$group,decreasing = T),]

x <- rt[,group$sample]
mat <- t(scale(t(x)))#归一化
df1 <- as.data.frame(mat)
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2

pdf(paste0('04.heatmap(',gse,').pdf'),  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("SKCM" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#104680", "white","#D80305"))(100))
dev.off()

png(paste0('04.heatmap(',gse,').png'),w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("SKCM" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#104680", "white","#D80305"))(100))
dev.off()

