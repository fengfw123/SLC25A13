rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "scRNA"
if(!dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)

if (! dir.exists("./00_rawdata/")){
  dir.create("./00_rawdata/")
}
setwd("./00_rawdata/")

library(pacman)
pacman::p_load(GEOquery, data.table, tidyverse, Seurat, patchwork,clusterProfiler, readr)
## PART A导入数据 --------------------------------------------------------------
folders <- list.files('./GSE222446_RAW/')
folders
dir <- paste0("GSE222446_RAW/", folders)
dir
scRNAlist <- list()
for(i in 1:length(dir)){
  raw_data<-  read.table(dir[i],header = T)#
  # colnames(raw_data)[1]<-'SYMBOL'
  counts <- raw_data
  counts <- counts[!duplicated(rownames(counts)),]
  # rownames(counts)<-counts$SYMBOL
  # counts <- counts[,c(2:ncol(counts))]
  # project <- tstrsplit(folders, "_", fixed = TRUE)[[1]]
  project <- str_split((str_split(dir[i],'[_]')[[1]][2]),'/')[[1]][2]   ##样本名称
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,
                                       min.cells = 3, min.features = 200,project=project)   ##将矩阵转换为Seurat对象，并对数据进行过滤
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = project)
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") }
  #计算核糖体基因比例
  # if(T){
  #   scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  # }
}
scRNA <- merge(x=scRNAlist[[1]],
               scRNAlist[2:length(scRNAlist)])
###线粒体基因表达占比(正则表达式，表示以MT-开头；scRNA[["percent.mt"]]这种写法会在meta.data矩阵加上一列)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")  
# ###计算核糖体基因比例
# scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

head(scRNA@meta.data)
system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))

#  PART B QC质控 ---------------------------------------------------------------
rm(list = ls())
setwd(('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/'))
if (! dir.exists("./01_QC")){
  dir.create("./01_QC")
}
setwd("./01_QC")
load('../00_rawdata/01.scRNA_orig.Rdata')


### 质控前 ###
paste0("质控前共有",length(colnames(scRNA)),"个细胞和", length(rownames(scRNA)),"个基因")   
# "质控前共有87449个细胞和26901个基因"

#绘图风格
theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  text = element_text(family = "Times"))

# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   


# 质控前小提琴图
violin <- VlnPlot(scRNA, 
                  pt.size = 0,
                  group.by = 'orig.ident',
                  features = plot.featrures)
# violin <- wrap_plots(plots = plots, nrow=1,ncol = 3,cols = palettes(category = "random",22,show_col=F))
violin


ggsave("01.vlnplot_before_qc.pdf", plot = violin, width = 10, height = 6)
ggsave("01.vlnplot_before_qc.png", plot = violin, width = 10, height = 6,dpi = 600)


### 质控后 ###
## 设置质控标准
minGene <- quantile(scRNA$nFeature_RNA,.02)   # 433
maxGene <- quantile(scRNA$nFeature_RNA,.98)   # 3059
maxUMI <- quantile(scRNA$nCount_RNA,.95)      # 10224


#  设置阈值
pctMT <- 20
minGene <- 400 #最小基因
maxGene <- 3000 # 最大基因
maxUMI <- 5000 # 最大细胞
minUMI <- 200 # 最小细胞

scRNA <- subset(scRNA, subset =
                  nCount_RNA > minUMI &
                  nCount_RNA < maxUMI &
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)
paste0("质控后共获得",length(colnames(scRNA)),"个细胞和", length(rownames(scRNA)), "个基因")   
# "质控后共获得33845个细胞和26901个基因"



# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
#GSE = "Patients_ID"

# 质控后小提琴图
# plots <- list()
# for(i in seq_along(plot.featrures)){
  # plots[[i]] = VlnPlot(scRNA, 
                       # cols = c("#ed556a"),
                       # pt.size = 0,
                       # features = plot.featrures[i]) + theme.set + NoLegend()}
# violin2 <- wrap_plots(plots = plots, nrow=1,ncol = 3)
violin2 <- VlnPlot(scRNA, 
                   pt.size = 0,
                   group.by = 'orig.ident',
                   features = plot.featrures)
violin2
dev.off()
ggsave("02.vlnplot_after_qc.pdf", plot = violin2, width = 10, height = 6)
ggsave("02.vlnplot_after_qc.png", plot = violin2, width = 10, height = 6,dpi = 600)
system.time(save(scRNA, file = "03.scRNA_qc.Rdata"))  


# PART C Normalize 标准化  -----------------------------------------------------
rm(list = ls())
setwd(('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/'))
if (! dir.exists("./02_normalization")){
  dir.create("./02_normalization")
}
setwd("./02_normalization")
load('../01_QC/03.scRNA_qc.Rdata')

#  合并count
scRNA <- JoinLayers(scRNA)

#标准化
scRNA.norm <- NormalizeData(scRNA,normalization.method = "LogNormalize",
                            scale.factor = 10000)

## 缩放数据
## 筛选2000个高变基因
all.genes <- rownames(scRNA.norm)
scRNA.norm <- FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA.norm), 10)


#  plot variable features with and without labels
#  保存基因特征方差图
pdf(file="01.feature_selection.pdf",width=10,height=6,family='Times',onefile=F)             
theme.set= theme(#axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 10,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))
plot1 <- VariableFeaturePlot(object = scRNA.norm)+theme.set
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
print(plot1+plot2)
dev.off()

#保存基因特征方差图
png(file="01.feature_selection.png",width=10,height=6,family='Times',units='in',res=600)              
plot1 <- VariableFeaturePlot(object = scRNA.norm)+theme.set
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 <- VariableFeaturePlot(object = scRNA.norm)+theme.set
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
print(plot1+plot2)
dev.off()

system.time(save(scRNA.norm, file = "02.scRNA_normalization.Rdata")) 


# PART D PCA  数据降维  --------------------------------------------------------
rm(list = ls())
setwd(('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/'))
if (!dir.exists("03_PCA")) {dir.create("03_PCA")}
setwd('03_PCA')

load('../02_normalization/02.scRNA_normalization.Rdata')
#   PCA降维之前的标准预处理步骤,这里对高可变基因进行数据标准化的操作,归一化表达量
scRNA.nor.sca <- ScaleData(scRNA.norm)

#   PCA分析   npcs = 20
scRNA.norm.pca <- RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca),npcs = 50 )
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
##
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
DimPlot(scRNA.norm.pca, reduction = "pca")
DimHeatmap(scRNA.norm.pca, dims = 1:5, cells = 2000, balanced = TRUE)


theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  text = element_text(family = "Times"))

##ElbowPlot,碎石图
pdf(file="01.ElbowPlot.pdf",width=8,height=6)
ElbowPlot(scRNA.norm.pca,reduction="pca",ndims = 50) +theme.set #根据ElbowPlot确定PC数量
dev.off()

png(file="01.ElbowPlot.png",width=8,height=6,units='in',res=600)
ElbowPlot(scRNA.norm.pca,reduction="pca",ndims = 50)+theme.set  
dev.off()


#  每个PC的p值分布和均匀分布,选择PC的数目，并且分别都对应不同的可视化
#  首先找到最佳聚类数
#  定义数据集的维度
#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50) #运行速度慢，影响JackStrawPlot
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)


##根据下面两张图可以选择PC数量
#JackStrawPlot
pdf(file="02.JackStrawPlot.pdf",width=12,height=6)
JackStrawPlot(object = scRNA.norm.pca, dims = 1:50,reduction = "pca")+theme.set
dev.off()

png(file="02.JackStrawPlot.png",width=12,height=6,units='in',res=600)
JackStrawPlot(object = scRNA.norm.pca, dims = 1:50,reduction = "pca")+theme.set
dev.off()


##主成分分析图形
pdf(file="03.PCA.pdf",width=8,height=6)
# DimPlot(object = scRNA.norm.pca, reduction = "pca",cols = c("DarkOrange", "SeaGreen"))+theme.set
DimPlot(object = scRNA.norm.pca, reduction = "pca")+theme.set
dev.off()

png(file="03.PCA.png",width=8,height=6,units='in',res=600)
# DimPlot(object = scRNA.norm.pca, reduction = "pca",cols = c("DarkOrange", "SeaGreen"))+theme.set
DimPlot(object = scRNA.norm.pca, reduction = "pca")+theme.set
dev.off()

# setwd(('/data/nas1/yuanshilin/project/03.YQWLMQ-10205-1/21_scRNA/03_PCA/'))
system.time(save(scRNA.norm.pca, file = "scRNA_pca.Rdata")) 

# PART E UMAP降维    -----------------------------------------------------------
rm(list = ls())
setwd(('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/'))
if (!dir.exists("04_UMAP")) {dir.create("04_UMAP")}
setwd('04_UMAP')
load('../03_PCA/scRNA_pca.Rdata')

pacman::p_load(tidyverse, scales, RColorBrewer, ggsci, wesanderson, IOBR)


pcSelect <- 30
scRNA.norm.umap <- FindNeighbors(object = scRNA.norm.pca, dims = 1:pcSelect)  #计算邻接距离
##对细胞分组,优化标准模块化
##这里的0.5决定了细胞亚群的多少，Seurat官网推荐0.4-1.2都可以进行尝试。数字越小，细胞亚群越少

##resolution参数表达聚类的分辨率，这个值大于0，一般都是在0-1范围里面调整，越大得到的cluster越多
#resolution可以反复调整，并不会改变降维的结果
scRNA.norm.umap <- FindClusters(object = scRNA.norm.umap, resolution = 0.5)                  

##UMAP聚类
scRNA.norm.umap  <- RunUMAP(scRNA.norm.umap , dims = 1:pcSelect)

##所有样本可视化
theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 20,face='bold',family = "Times"),
  text = element_text(family = "Times"))

pdf(file="01.UMAP_all.pdf",width=10,height=8,family='Times')
DimPlot(scRNA.norm.umap ,
        reduction = "umap",
        label.size = 6,label = TRUE,pt.size=1.0, 
        cols = palettes(category = "random",29,show_col=F))+theme.set
dev.off()

png(file="01.UMAP_all.png",width=10,height=8,family='Times',units='in',res=600)
DimPlot(scRNA.norm.umap ,
        reduction = "umap",
        label.size = 6,label = TRUE, pt.size=1.0,
        cols = palettes(category = "random",24,show_col=F))+theme.set
dev.off()

# pdf(file="02.UMAP_Group.pdf",width=14,height=8,family='Times')
# DimPlot(scRNA.norm.umap ,
#         reduction = "umap",
#         label.size = 6,label = T,pt.size=1.0,
#         split.by  = "group",
#         cols = palettes(category = "random",24,show_col=F))+theme.set
# dev.off()
# 
# png(file="02.UMAP_Group.png",width=14,height=8,family='Times',units='in',res=600)
# DimPlot(scRNA.norm.umap ,
#         reduction = "umap",
#         label.size = 6,label = T, pt.size=1.0,
#         split.by  = "group",
#         cols = palettes(category = "random",22,show_col=F))+theme.set
# dev.off()

system.time(save(scRNA.norm.umap, file = "scRNA_UMAP.Rdata")) 


# PART F 细胞注释    -----------------------------------------------------------
rm(list = ls())
setwd('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/')
if (! dir.exists("./05_Annotation")){
  dir.create("./05_Annotation")
}
setwd("./05_Annotation")

load('../04_UMAP/scRNA_UMAP.Rdata')
pacman::p_load(SingleR, celldex, reshape2, IOBR,tidyverse, colorspace)

               

scRNA_nor_singleR <- GetAssayData(scRNA.norm.umap, layer = "data")

hpca.se <- HumanPrimaryCellAtlasData()
#hpca.se<-BlueprintEncodeData()
#hpca.se<-DatabaseImmuneCellExpressionData()
clusters <- scRNA.norm.umap@meta.data$seurat_clusters
pred.hesc <- SingleR(test = scRNA_nor_singleR,
                     ref = hpca.se,
                     labels = hpca.se$label.main,
                     method = "cluster", 
                     clusters = clusters, 
                     quantile = 0.8,
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")

table(pred.hesc$labels)
### 9
celltype <- data.frame(ClusterID=rownames(pred.hesc), 
                       celltype=pred.hesc$labels, 
                       stringsAsFactors = F) 

## 手动注释
DefaultAssay(scRNA.norm.umap) <- "RNA"

### 不能有重复
list_gene <- list(
  Melanoma = c('MLANA','ITGAM','PRAME','SOX10','S100A1'),
  Tcell = c('CD2',	'CD3D',	'CD3E', 'CD4', 'CD8A', 'CD8B'),
  Treg = c('FOXP3', 'IL2RA'),
  NK = c('GNLY', 'KLRD1','NKG7'),
  Bcell = c('CD79A','MS4A1', 'BANK1'),
  PlasmaB = c('JCHAIN','MZB1', 'IGHG1', 'IGHG4'),
  pDC = c('LILRA4','LAMP5','PLD4'),
  Myeloid = c('CD14','CD68','CSF1R','MRC1','LYZ','CPVL','G0S2'),
  Fibroblast = c('EPCAM', 'KRT19','DCN', 'FGF7'),
  # MC = c('ACTA2', 'MYH11','TAGLN','CNN1'),
  Endothelial = c("PECAM1", "VWF", "CD34", "CDH5")
)


theme.set = theme(
  axis.title = element_text(size = 18, face = "bold", family = "Times"),
  axis.text.x = element_text(angle = 45, vjust = 1,size = 14, family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))


# experssion dot plot 
Idents(scRNA.norm.umap) <- "seurat_clusters"
DotPlot(scRNA.norm.umap, features = list_gene) + RotatedAxis() + theme.set +
  scale_color_gradientn(colours = c('#0f59a4','white','#d11a2d'))

pdf(file = "01.dotplot.pdf", width = 24, height = 12, family = 'Times')
DotPlot(scRNA.norm.umap, features = list_gene) + RotatedAxis() + theme.set+
  scale_color_gradientn(colours = c('#0f59a4','white','#d11a2d'))
dev.off()

png(file = "01.dotplot.png", width = 24, height = 12, family = 'Times', units = 'in', res = 600)
DotPlot(scRNA.norm.umap, features = list_gene) + RotatedAxis() + theme.set+
  scale_color_gradientn(colours = c('#0f59a4','white','#d11a2d'))
dev.off()

### 手动给cluster分组
cluster.ids = c(
  "0" = "Melanoma",  
  "1" = "Tcell",  
  "2" = "Tcell",  
  "3" = "Bcell", 
  "4" = "Melanoma",   
  "5" = "Tcell",
  "6" = "Melanoma",
  "7" = "Tcell",
  "8" = "Myeloid",
  "9" = "Tcell",
  "10" = "Melanoma",
  "11" = "Melanoma",
  "12" = "Tcell",
  "13" = "Tcell",
  "14" = "Melanoma",
  "15" = "pDC",
  "16" = "Myeloid",
  "17" = "Endothelial",
  "18" = "Melanoma",
  "19" = "Fibroblast",
  "20" = "Melanoma",
  "21" = "Tcell"
)

scRNA.norm.umap <- RenameIdents(scRNA.norm.umap, cluster.ids)
scRNA.norm.umap$celltype <- Idents(scRNA.norm.umap)
table(scRNA.norm.umap@meta.data$celltype)

### 绘制注释umap图
theme.set2 = theme(
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12, face = "bold", family = "Times"),
  axis.text.y = element_text(size = 12, face = "bold", family = "Times"),
  legend.text = element_text(size = 14, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

pdf(file="02.celltype.pdf",width=9,height=6,family='Times')
DimPlot(scRNA.norm.umap ,reduction = "umap",
        label.size = 6,label.color = 'black',label = T,
        cols = palettes(category = "random",20,show_col=F))+theme.set2
dev.off()

png(file="02.celltype.png",width=9,height=6,family='Times',units='in',res=600)
DimPlot(scRNA.norm.umap ,reduction = "umap",
        label.size = 6,label.color = 'black',label =T,
        cols = palettes(category = "random",20,show_col=F))+theme.set2
dev.off()


###############################我是一条分割线################################### 另一种方法
#  可以在聚类图上查找
# ###  Macrophage ('CD68', 'ITGAM', 'ITGAX', 'CD14', 'CSF1R'),
# p1 <-  FeaturePlot(scRNA.norm.umap,
#                  features = c('CD68', 'ITGAM', 'ITGAX', 'CD14', 'CSF1R'),
#                  repel = T, label = T, slot = "data", 
#                  cols = c("lightgrey","red"), ncol = 3, order = T)
# p1
# # ###  myeloid/macrophage (CD14, CSF1R and AIF1),
# p2  <-  FeaturePlot(scRNA.norm.umap, 
#                     features = c("LILRA4","CLEC4C",'GZMB'),
#                     repel = T, label = T, slot = "data",
#                     cols = c("lightgrey","red"), ncol = 3, order = T)
# p2
##11
###############################我是一条分割线################################### 另一种方法

## marker基因在细胞中表达情况

all.markers <- as.data.frame(t(sapply(list_gene, "[", i = 1:max(sapply(list_gene, length)))))

write.csv(all.markers,file = '03.AllMarkers.csv',row.names = F,quote = F)

marker_gene <- unlist(list_gene)

pdf(file = '04.marker.exp.pdf',width=24,height=12,family='Times')
DotPlot(scRNA.norm.umap,features = list_gene)+
  RotatedAxis()+theme.set+
  scale_color_gradientn(colours = c('#0f59a4','#e9ccd3','#d11a2d'))
dev.off()
png(file = '04.marker.exp.png',width=24,height=12, units = 'in',family='Times', res = 600)
DotPlot(scRNA.norm.umap,features = list_gene)+
  RotatedAxis()+theme.set+
  scale_color_gradientn(colours = c('#0f59a4','#e9ccd3','#d11a2d'))
dev.off()


system.time(save(scRNA.norm.umap, file = "scRNA.umapcell.type.Rdata"))

##换个形式画一下比例图
all.count <- data.frame(celltype=scRNA.norm.umap$celltype)
table(all.count$celltype)
# 10

plot.celltype <- table(all.count$celltype) %>% as.data.frame() %>% recast(Var1)
colnames(plot.celltype)[1] <- "Cell"
plot.celltype <- cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
colnames(plot.celltype)[2:3] <- c("Count_SKCM","Proportion_SKCM")
plot.celltype <- plot.celltype[order(plot.celltype$Count_SKCM,decreasing = T),]
write.csv(plot.celltype, "05.CellCount.csv", row.names = F, quote = F)

dat.celltype <- plot.celltype[c(1,3)]

colnames(dat.celltype)[2] = c("SKCM")
dat.celltype[2] <- 100 * dat.celltype[2]

dat.plot <- melt(dat.celltype, id.vars = "Cell", variable.name = "Sample", value.name = "value")
dat.plot$Cell <- factor(dat.plot$Cell, levels = plot.celltype$Cell)

dat.plot$label = round(dat.plot$value, 2) %>% paste0(.,"%")
ggplot(data = dat.plot, mapping = aes(x = value, y = Cell)) + 
  geom_text(aes(x = value + 1, label = label), hjust = 0) +
  geom_bar(aes(fill = Cell), stat = "identity") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values = palettes(category = "random",22,show_col=F)) + 
  xlim(c(0,100)) +
  xlab("Cell Fraction") + ylab("Cell Type") +
  theme_minimal(base_size = 16) + 
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"),
        legend.position = "none", 
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, colour = 'black',family = "Times"),
        axis.text.y = element_text(size = 12, colour = 'black',family = "Times"),
        axis.title = element_text(size = 16, face = "bold", family = "Times"))

ggsave("06.CellProportion.pdf", width = 4, height = 4)
ggsave("06.CellProportion.png", width = 4, height = 4, units = "in", dpi = 600)


# PART G 关键细胞    -----------------------------------------------------------
rm(list = ls())
setwd(('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/'))
if (! dir.exists("./06_hubcell")){dir.create("./06_hubcell")}
setwd("./06_hubcell")

###---- load data ----###
load('../05_Annotation/scRNA.umapcell.type.Rdata')

# load gene file
gene <- read.csv(file = "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/06_Model101/Results/03.fea_df.csv")
gene <- data.frame(symbol = gene[gene$algorithm == "StepCox[backward]+RSF",]$features)

# add module score find hub cell
# sce.all <- AddModuleScore(object = scRNA.norm.umap, list(gene$symbol))
# head(sce.all@meta.data)

# experssion list 
expr.list <- list()

# extract experssion data 
for (i in c(1:nrow(gene))) {
  vln.df <- data.frame(exp = scRNA.norm.umap@assays$RNA$data[gene[i,],])
  vln.df$sample <- rownames(vln.df)
  vln.df <- reshape::melt(vln.df, id = "sample")
  vln.df <- filter_if(vln.df, is.numeric, all_vars((.)>0))
  anno <- data.frame(cell = scRNA.norm.umap@meta.data[,c("celltype")],
                     sample = rownames(scRNA.norm.umap@meta.data))
  vln.df <- inner_join(vln.df, anno, by = "sample")
  vln.df$symbol <- paste0(gene[i,])
  expr.list[[i]] <- vln.df
}

# experssion data 
expr.data <- do.call(rbind, expr.list) %>% dplyr::select("symbol", "cell", "value")

# bar fig 
bar.plot <- ggplot(expr.data, aes(x = symbol, y = value, fill = cell)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gene Expression by Cell Type",
       x = "Gene",
       y = "Expression Level") +
  theme_classic() +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494"))+
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", colour = 'black',family = "Times"),
        axis.text.x = element_text(size = 12, colour = 'black',family = "Times"),
        axis.text.y = element_text(size = 12, colour = 'black',family = "Times"),
        axis.title = element_text(size = 16, face = "bold", family = "Times"),
        # panel.border = element_blank(), # 去掉边框
        axis.line = element_line(color = "black"),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank())
ggsave("01.geneexpr_bar.pdf", bar.plot, width = 10, height = 8)
ggsave("01.geneexpr_bar.png", bar.plot, width = 10, height = 8, units = "in", bg = "white")


###---- choose hubcell ----####
# extract hubcell seurat 
scRNA.norm.hub <- subset(scRNA.norm.umap, idents = c("Melanoma","Tcell", "Bcell", "Myeloid"))

# hubgene
hubgene <- gene

# umap expr
p1 <- FeaturePlot(scRNA.norm.hub,
                  features = hubgene$symbol, repel = T, label = T,
                  slot = "data", cols = c("lightgrey","#D80305"), ncol = 2, order = T)
p1
ggsave(filename = '02.cell_map.pdf', p1, width = 10, height = 7, units = "in")
ggsave(filename = '02.cell_map.png', p1, width = 10, height = 7, units = "in", dpi = 600, bg = "white")


###---- pseudotime anlysis ----####
rm(list = ls())
setwd('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/')
dir <- "07_pseudotime"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

###---- load packages ----###
library(tidyverse); library(scales); library(RColorBrewer); library(ggsci); library(wesanderson)
library(IOBR); library(Seurat)
devtools::load_all('/data/nas1/yuanshilin/Database/scRNA/monocle/')
library(monocle)

###---- lpseudotime ----###
# load seurat data 
load("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/05_Annotation/scRNA.umapcell.type.Rdata")

# hucells
hub.cell <- c("Melanoma","Tcell", "Bcell", "Myeloid")

# load gene file
hubgene <- read.csv(file = "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/06_Model101/Results/03.fea_df.csv")
hubgene <- data.frame(symbol = hubgene[hubgene$algorithm == "StepCox[backward]+RSF",]$features)

# nogroup cells
nogroup <- c("Melanoma","Tcell", "Bcell", "Myeloid")

i <- 1
for (i in c(1:length(hub.cell))) {
  setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/07_pseudotime/")
  tmp.dir <- hub.cell[i]
  if (! dir.exists(tmp.dir)){dir.create(tmp.dir)}
  setwd(tmp.dir)
  # dimensionality reduction
  scRNA.single.cell <- scRNA.norm.umap[ ,scRNA.norm.umap$celltype %in% c(hub.cell[i])]
  scRNA.single.cell <- NormalizeData(scRNA.single.cell,
                                     normalization.method = "LogNormalize", scale.factor = 10000)
  scRNA.single.cell <- FindVariableFeatures(scRNA.single.cell,
                                            selection.method = 'vst', nfeatures = 1000)
  scRNA.single.cell <- ScaleData(scRNA.single.cell)
  scRNA.single.cell <- RunPCA(scRNA.single.cell, features = VariableFeatures(object = scRNA.single.cell))
  scRNA.single.cell <- FindNeighbors(scRNA.single.cell, dims = 1:20)
  scRNA.single.cell <- FindClusters(scRNA.single.cell, resolution = 0.4)
  scRNA.single.cell <- RunUMAP(scRNA.single.cell, dims = 1:20)
  
  theme.set <- theme(
    plot.title = element_text(size = 20,face = "bold", family = "Times"),
    axis.title = element_text(size = 18, face = "bold", family = "Times"),
    axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
    axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
    legend.text = element_text(size = 16, face = "bold", family = "Times"),
    text = element_text(family = "Times"))
  p <- DimPlot(scRNA.single.cell ,
               reduction = "umap",
               group.by = "seurat_clusters",
               label.size = 6, label = TRUE, pt.size=1.0,
               cols = palettes(category = "random", 14, show_col = F)) + theme.set
  ggtitle(paste0('Subtype of ', hub.cell[i], ' cells'))
  ggsave(filename = paste0('01.subtype_,', hub.cell[i], '.pdf'), p, width = 6, height = 5)
  ggsave(filename = paste0('01.subtype_,', hub.cell[i], '.png'), p, width = 6, height = 5, dpi = 600, bg = "white")
  system.time(save(scRNA.single.cell, file = paste0("scRNA.", hub.cell[i], ".Rdata")))
  
  # pseudotime anlysis
  single.matrix <- GetAssayData(scRNA.single.cell, layer = "count", assay = "RNA")
  feature.ann <- data.frame(gene_id = rownames(single.matrix), gene_short_name = rownames(single.matrix))
  
  # gene feature
  rownames(feature.ann) <- rownames(single.matrix)
  
  # AnnotatedDataFrame
  single.fdat <- new("AnnotatedDataFrame", data = feature.ann)
  
  sample.ann <- scRNA.single.cell@meta.data
  cell.type <- Idents(scRNA.single.cell) %>% data.frame(.)
  sample.ann <- cbind(sample.ann, cell.type)
  
  # add subcluster cell information
  colnames(sample.ann)[ncol(sample.ann)] = "Subtype"
  
  single.pdat <- new("AnnotatedDataFrame", data = sample.ann)
  
  single.cds <- newCellDataSet(single.matrix, phenoData = single.pdat, featureData = single.fdat,
                               expressionFamily = negbinomial.size())
  single.cds <- estimateSizeFactors(single.cds)
  single.cds <- estimateDispersions(single.cds)
  
  # Free resources that aren't needed 
  WGCNA::collectGarbage()
  
  # add experssion
  express.genes <- VariableFeatures(scRNA.single.cell)
  single.exp.cds <- setOrderingFilter(single.cds, express.genes)
  
  single.exp.cds <- reduceDimension(single.exp.cds, max_components = 2,
                                    verbose = T, norm_method = "log")   #, residualModelFormulaStr = "~nFeature_RNA+group"
  single.exp.cds <- orderCells(single.exp.cds)
  save(single.exp.cds, file = paste0(hub.cell[i],"_cds_order.Rdata"))
  
  p1 <- plot_cell_trajectory(single.exp.cds, color_by = "Pseudotime", cell_size = 1, theta = 180,
                             size = 1, show_backbone=TRUE, show_branch_points = F) +
    theme(text = element_text(size = 14))
  p1
  ggsave(paste0("02.", hub.cell[i],"_Trajectory.png"), p1, width = 6, height = 6, units = "in", dpi = 600)
  ggsave(paste0("02.", hub.cell[i],"_Trajectory.pdf"), p1, width = 6, height = 6, units = "in")
  
  
  p2 <- plot_cell_trajectory(single.exp.cds, color_by = "seurat_clusters", cell_size = 1, theta = 180,
                             size = 1, show_backbone=TRUE, show_branch_points = F) +
    scale_color_simpsons()+
    theme(text = element_text(size = 14))
  p2
  ggsave(paste0("03.", hub.cell[i], "_type_Trajectory.png"), p2, width = 6, height = 6, units = "in", dpi = 600)
  ggsave(paste0("03.", hub.cell[i], "_type_Trajectory.pdf"), p2, width = 6, height = 6, units = "in")
  
  p3 <- plot_cell_trajectory(single.exp.cds, color_by = "State", cell_size = 1, theta = 180,
                             size = 1, show_backbone = TRUE, show_branch_points = F) +
    scale_color_npg() +
    theme(text = element_text(size = 14))
  p3
  ggsave(paste0("04.", hub.cell[i], "_State_Trajectory.png"), p3, width = 6, height = 6, units = "in", dpi = 600)
  ggsave(paste0("04.", hub.cell[i], "_State_Trajectory.pdf"), p3, width = 6, height = 6, units = "in")
  
  if(!hub.cell[i] %in% nogroup){
    p4 <- plot_cell_trajectory(single.exp.cds, color_by = "seurat_clusters", cell_size = 1, theta = 180,
                               size = 1, show_backbone = TRUE, show_branch_points = F) +
      scale_color_simpsons()+
      theme(text = element_text(size = 14))+
      facet_wrap(~ Group, nrow = 1)
    p4
    ggsave(filename = '05.celltype_Trajectory(group).pdf', p4, width = 10, height = 6, units = "in", dpi = 600)
    ggsave(filename = '05.celltype_Trajectory(group).png', p4, width = 10, height = 6)
    
    # cell experssion for pseudotime plot
    cds.subset <- single.exp.cds[hubgene$symbol ,]
    theme.set <-   theme(legend.position = "right")+
      theme(axis.title.x = element_text(size = 22, color = 'black', face = "bold", family = 'Times'),
            axis.text.x = element_text(size = 18, family = 'Times'),
            axis.title.y = element_text(size = 22,color='black', face = "bold", family = 'Times'),
            axis.text.y = element_text(size = 18,  family = 'Times'),
            legend.title=  element_text(size = 20, color = 'black', face = "bold", family = 'Times'),
            legend.text = element_text(size = 18,  face = "bold", family = 'Times'),
            title = element_text(size = 20, color = 'black', face = "bold", family = 'Times'))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    p1 <- plot_genes_in_pseudotime(cds.subset , color_by = "Pseudotime") + theme.set
    
    p2 <- plot_genes_in_pseudotime(cds.subset , color_by = "State") + theme.set
    
    p3 <-plot_genes_in_pseudotime(cds.subset , color_by = "Group") + theme.set
    plot <- p1|p2|p3
    
    ggsave(paste0('06.trajectory_gene_', hub.cell[i], '.pdf'), plot, height = 8, width = 18, family='Times')
    ggsave(paste0('06.trajectory_gene_', hub.cell[i], '.png'), plot, height = 8, width = 18, dpi = 600)
  }else{
    cds.subset <- single.exp.cds[hubgene$symbol ,]
    theme.set <-   theme(legend.position = "right")+
      theme(axis.title.x = element_text(size = 22, color = 'black', face = "bold", family = 'Times'),
            axis.text.x = element_text(size = 18, family = 'Times'),
            axis.title.y = element_text(size = 22,color='black', face = "bold", family = 'Times'),
            axis.text.y = element_text(size = 18,  family = 'Times'),
            legend.title=  element_text(size = 20, color = 'black', face = "bold", family = 'Times'),
            legend.text = element_text(size = 18,  face = "bold", family = 'Times'),
            title = element_text(size = 20, color = 'black', face = "bold", family = 'Times'))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    p1 <- plot_genes_in_pseudotime(cds.subset , color_by = "Pseudotime") + theme.set
    
    p2 <- plot_genes_in_pseudotime(cds.subset , color_by = "State") + theme.set
    
    plot <- p1|p2
    
    ggsave(paste0('05.trajectory_gene_', hub.cell[i], '.pdf'), plot, height = 8, width = 12, family='Times')
    ggsave(paste0('05.trajectory_gene_', hub.cell[i], '.png'), plot, height = 8, width = 12, dpi = 600)
  }
}


####---- cell chat ----####
rm(list = ls())
setwd('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/')
dir <- "08_cell_chat"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

###---- load packages ----###
library(CellChat); library(Seurat); library(ggplot2); library(tidyverse)

####---- load data ----###
# load seurat data 
load("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/13_scRNA/05_Annotation/scRNA.umapcell.type.Rdata")

# hubcell seurat
scRNA.hub.cell <- subset(scRNA.norm.umap, idents = c("Melanoma","Tcell", "Bcell", "Myeloid"))

# normalized data matrix
data.input <- scRNA.hub.cell@assays$RNA$data 
meta <- scRNA.hub.cell@meta.data # a dataframe with rownames containing cell mata data
table(meta$Group)

###---- disease cell ----###
# data set
di.dat <- data.input
di.meta <- meta

# CellChat need matrix data and meta 
di.cellchat <- createCellChat(object = di.dat, meta = di.meta, group.by = "celltype")

# add meta data
di.cellchat <- addMeta(di.cellchat, meta = di.meta)

# change labels as your needed
di.cellchat <- setIdent(di.cellchat, ident.use = "celltype") 
levels(di.cellchat@idents)
unique(di.cellchat@idents)

di.cellchat@idents <- droplevels(di.cellchat@idents,
                                 exclude = setdiff(levels(di.cellchat@idents),
                                                   unique(di.cellchat@idents)))

# cell numbers
group.size <- as.numeric(table(di.cellchat@idents)) 

# human or mouse
CellChatDB <- CellChatDB.human
# CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.interaction <- CellChatDB$interaction

# use “Secreted Signaling” anylsis 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
di.cellchat@DB <- CellChatDB.use

# The expression data is further preprocessed to save computing power
di.cellchat <- subsetData(di.cellchat)

# overexpression genes (ligand-receptor) are identified first
di.cellchat <- identifyOverExpressedGenes(di.cellchat)

# the overexpressed interactions between the overexpressed ligands and receptors are then identified
di.cellchat <- identifyOverExpressedInteractions(di.cellchat)

# project gene expression data onto PPI network (optional) 
# tip: human or mouse 
# cellchat <- projectData(cellchat, PPI.human)
di.cellchat <- projectData(di.cellchat, PPI.mouse)

# Calculate the probability of intercellular communication and predict the communication network
di.cellchat <- computeCommunProb(di.cellchat)

# output net 
di.net <- subsetCommunication(di.cellchat) 
write.csv(di.net,file = '01.net.disease.csv', row.names = F, quote = F)

# the level of signaling pathways further infers intercellular communication and calculates aggregation networks
di.cellchat <- computeCommunProbPathway(di.cellchat)
di.cellchat <- aggregateNet(di.cellchat)
group.size <- as.numeric(table(di.cellchat@idents))

# cell chat heatmap
netVisual_heatmap(di.cellchat)
pdf('02.number_of_interactions_disease.pdf', width = 6, height = 5, family = 'serif')
netVisual_heatmap(di.cellchat)
dev.off()
png('02.number_of_interactions_disease.png', width = 500, height = 400, family = 'serif')
netVisual_heatmap(di.cellchat)
dev.off()

# save rds 
saveRDS(di.cellchat, file = 'diseasecellchat.rds')

# cell chat net of count 
pdf(file = '03.net_number_disease.pdf', width = 7, height = 9, family = 'serif')
netVisual_circle(di.cellchat@net$count,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()

png(filename = '03.net_number_disease.png', width = 7, height = 9, units = 'in', res = 600, family = 'serif')
netVisual_circle(di.cellchat@net$count,
                 vertex.weight = group.size, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
dev.off()

# cell chat net of weight
pdf(file = '04.net_weight_disease.pdf', width = 7, height = 9, family = 'serif')
netVisual_circle(di.cellchat@net$weight,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()

png(filename = '04.net_weight_disease.png',w=7,h=9,units='in',res=600, family='serif')
netVisual_circle(di.cellchat@net$weight,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()

levels(di.cellchat@idents)
# show all the significant interactions (L-R pairs)
# dot plot
di.cellchat@data.signaling
pdf(file = '05.buble_disease.pdf', width = 14, height = 9)
netVisual_bubble(di.cellchat, remove.isolate = FALSE)
dev.off()
png('05.buble_disease.png', width = 14, height = 9, units = 'in', res = 600)
netVisual_bubble(di.cellchat, remove.isolate = FALSE)
dev.off()








