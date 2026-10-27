rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "11_GSEA"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)

###--- load packages ----####
library(data.table); library(org.Hs.eg.db); library(clusterProfiler); library(biomaRt)
library(enrichplot); library(limma); library(tidyverse); library(lance);library(GSVA)
library(GSEABase); library(DESeq2); library(dplyr)

###---- load data ----####
# count data
dat.tcga <- read.csv("../00_rawdata/SKCM/02.dat.tcga.csv", row.names = 1)
colnames(dat.tcga) <- gsub('.','-',colnames(dat.tcga), fixed = T)

# load group
group <- read.csv("../07_Riskmodel/01.train_riskScore.csv", row.names = 1) %>% dplyr::select(c("id","risk"))
group <- group[order(group$risk),]
gse.exp <- dat.tcga[,group$id]

# dds data and group 
colData <- group
colnames(colData) <- c('sample','group')
colData$group <- factor(colData$group, levels = c('High','Low'))
gse.exp <- round(gse.exp, 0)

# different anlysis
dds <- DESeqDataSetFromMatrix(countData = gse.exp, colData = colData, design = ~group)
dds <- dds[rownames(counts(dds)) > 1,]
dds <- DESeq(dds)
dds

# different result
res <- results(dds, contrast = c("group","High","Low"))
res <- res[order(res$padj), ]
head(res)
summary(res)
table(res$padj < 0.05)
allGeneSets <- as.data.frame(res) %>%  na.omit()

logFCcutoff <- 1
allGeneSets$change = as.factor(
  ifelse(allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff,
         ifelse(allGeneSets$log2FoldChange > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > 1)
genelist <- allGeneSets$log2FoldChange
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
dim(DEGeneSets)
write.csv(allGeneSets,'01.allGeneSet.csv',row.names = T,quote = F)
write.csv(DEGeneSets,'02.DEGeneSet.csv',row.names = T,quote = F)

## GSEA KEGG-----
kegg_set<- read.gmt("/data/nas1/yuanshilin/Database/Enrichment/c2.kegg.v7.4.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
write.csv(kegg_result,file = '03.GSEA_KEGG_result.csv',row.names = F,quote = F)

dim(kegg_result)
#trace(gseaplot2,edit=T)   #theme(legend.position = c(0.7, 0.8) 26行可以更改 前一个数值为左右 后一个数值为上下

p <- gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
               title = 'KEGG GSEA',
               base_size = 10,
               rel_heights = c(1, 0.2, 0.2))
p[[1]] <- p[[1]] + theme(legend.position = "right")          

pdf('04.GSEA(KEGG).pdf', width = 12, height = 6)
print(p)
dev.off()

png('04.GSEA(KEGG).png', width = 12, height = 6, units = "in", res = 600)
print(p)
dev.off()
