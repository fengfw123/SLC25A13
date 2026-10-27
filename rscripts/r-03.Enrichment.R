rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "03_Enrichment"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)

pacman::p_load(clusterProfiler, ggplot2,org.Hs.eg.db, tidyverse,GOplot,enrichplot,
               RColorBrewer, ggrepel)
################################ 富集分析 ######################################
# 01.原始数据 ------------------------------------------------------------------
gene.inter <- read.csv("../02_Hubgene/02.hubgenes.csv", row.names = 1)

database <- 'org.Hs.eg.db'
### 获得gene的Ensembol ID
gene_transform <- bitr(gene.inter$symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = database)

# 02.GO富集---------------------------------------------------------------------
### CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，
### 选自己需要的,一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = database,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                # qvalueCutoff = 0.05,
                readable = TRUE)
save(ego, file = 'GO.Rdata')
### go barplot -----------------------------------------------------------------
go <- ego # 留存结果
enrich_go <- mutate(go, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
go_result <- enrich_go@result
if(median(go_result$richFactor) < 0.1){go_result$richFactor <- go_result$richFactor*10}

### count排序 
GO_ALL <- go_result    # 防止后面出错
GO_ALL <- GO_ALL[GO_ALL$pvalue < 0.05,]
BP <- GO_ALL[GO_ALL$ONTOLOGY=='BP', ]
BP <- BP[order(BP$Count,decreasing = T),]
CC <- GO_ALL[GO_ALL$ONTOLOGY=='CC', ]
CC <- CC[order(CC$Count,decreasing = T),]
MF <- GO_ALL[GO_ALL$ONTOLOGY=='MF', ]
MF <- MF[order(MF$Count,decreasing = T),]

paste0("共得到",dim(GO_ALL)[[1]],"个GO富集结果，其中",dim(BP)[[1]],"个生物学过程，",dim(CC)[[1]],"个细胞组分，",dim(MF)[[1]],"个分子功能")
# "共得到 973 个GO富集结果，其中 770 个生物学过程，30 个细胞组分，173 个分子功能"
write.table(GO_ALL,file = "00.GO_result.csv",sep = ",",quote = F,row.names = F)
write.table(as.data.frame(BP), '00.GO_BP.csv', sep = ',', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(CC), '00.GO_CC.csv', sep = ',', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(MF), '00.GO_MF.csv', sep = ',', row.names = FALSE, quote = FALSE)

display_number = c(5, 5, 5)
ego_result_BP <- as.data.frame(BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(MF)[1:display_number[3], ]
paste0(ego_result_BP$Description,collapse = "（）；")
paste0(ego_result_CC$Description,collapse = "（）；")
paste0(ego_result_MF$Description,collapse = "（）；")
col1 <- c("#E64B35B2", "#DC0000B2",'#E63863','#B2182B', '#D6604D','#EB4B17')
col2 <- c("#4DBBD5B2", "#0073C2B2",'#008280B2','#4393C3', '#23452F','#68A180')
col3 <- c("#E4C755", "#9FA3A8",'#D6E7A3','#F4A371FF','#F9C187FF','#FACF94FF')
Col1 <- sample(col1,1)
Col2 <- sample(col2,1)
Col3 <- sample(col3,1)

if(T){
  go_enrich_df <- data.frame(
    ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
    Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
    GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
    type=factor(c(rep("BP:biological process", display_number[1]), 
                  rep("CC:cellular component", display_number[2]),
                  rep("MF:molecular function", display_number[3])), 
                levels=c("BP:biological process", "CC:cellular component","MF:molecular function" )))
  go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
  COLS <- c(Col1, Col2,Col3)
  
  p <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #横纵轴取值
    geom_bar(stat="identity", width=0.8) + #柱状图的宽度，可以自己设置
    scale_fill_manual(values = COLS) + ###颜色
    coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
    xlab("GO term") + 
    ylab("Gene Number") + 
    labs(title = "GO Terms")+
    theme_bw()+theme(
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
      legend.text = element_text(face="bold",color="black",family = "Times",size=12),
      plot.title = element_text(hjust = 0.5, face = "bold",color = "black",family = "Times",size = 18),
      axis.text.x = element_text(face = "bold",color = "black",size = 14),
      axis.text.y = element_text(face = "bold",color = "black",size = 14),
      axis.title.x = element_text(face = "bold",color = "black",family = "Times",size = 18),
      axis.title.y = element_text(face = "bold",color = "black",family = "Times",size = 18),
      plot.subtitle = element_text(hjust = 0.5,family = "Times", size = 14, face = "italic", colour = "black")
    )+scale_x_discrete(labels=function(x) str_wrap(x, width=40))
  ggsave(filename = "01.GO_barplot.pdf", height = 8, width = 12, p)
  ggsave(filename = "01.GO_barplot.png", height = 8, width = 12, p, dpi = 600)
}


### go circ --------------------------------------------------------------------
enrich_ID <- gene.inter
select_DERNA <- read.csv("../01_DEGs/02.DEGs_sig(GSE15605).csv", check.names=F, header=T, row.names = 1)

if(str_detect(colnames(select_DERNA)[1],'logFC')){
  colnames(select_DERNA)[1] <- 'log2FoldChange'
}
gene_symbol <- clusterProfiler::bitr(
  geneID = enrich_ID$x,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db")

go.circ = 1
if(go.circ == 1){
  go_tem <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
  go_dat <-  data.frame(Category = go_tem[,'ONTOLOGY'],ID = go_tem[,'ID'],Term = go_tem[,'Description'], 
                        Genes = gsub("/", ", ", go_tem[,'geneID']), adj_pval = go_tem[,'pvalue'])
  genelist <- data.frame(ID = rownames(select_DERNA), logFC = select_DERNA$log2FoldChange)
  row.names(genelist)=genelist[,1]
  library(GOplot)
  circ <- circle_dat(go_dat, genelist)
  go_tem$pvalue <- -log10(go_tem$pvalue)
  go_enrich_df <- data.frame(
    category = go_tem$ONTOLOGY, 
    gene_num.min = 0,
    gene_num.max = max(go_tem$Count),
    gene_num.rich = go_tem$Count,
    "-log10.p" = go_tem$pvalue,
    up.regulated = 0,
    down.regulated = 0,
    rich.factor = go_tem$richFactor,
    Description = go_tem$ID)
  circ$sig <- ifelse(circ$logFC > 0,"up","down")
  tem1 <- go_enrich_df$Description
  i =1
  for(i in c(1:length(tem1))){
    x = tem1[i]
    circ1 <- circ[grep(x,circ$ID),]
    num1 <- nrow(circ1)
    tem2 <- as.data.frame(table(circ1$sig))
    if(nrow(tem2) == 2){
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[2]]
    }else if(tem2[1,1] == "up"){
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[1]]
    }else{go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]}
  }
  
  go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
  rownames(go_enrich_df) <- go_enrich_df$Description
  pdf(file = paste0("02.GO_circplot.pdf"),width = 12,height = 6)
  a <- dev.cur()   
  png(file = paste0("02.GO_circplot.png"),width= 12, height= 6, units="in", res=300)
  dev.control("enable")
  par(mar = c(2,2,2,2),cex=1,family="Times")
  circle_size = unit(1, 'snpc')
  #加载 circlize 包
  library(circlize)
  ##整体布局
  circos.par(gap.degree = 2, start.degree = 90)
  ##第一圈，绘制 ko
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep(Col1, 5), rep(Col2, 5), rep(Col3, 5))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  # ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  # sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.4, major.tick.percentage = 0.4, labels.niceFacing = FALSE)  # 绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  # 将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(pvalue)")
  ##第二圈，绘制富集的基因和富集 p 值
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(pvalue)')]
  label_data <- go_enrich_df['gene_num.rich'] 
  p_max <- round(max(go_enrich_df$`-log10(pvalue)`)) + 1 
  colorsChoice <- colorRampPalette(c('white','#4393C3'))
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.3, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  
  
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('#FDDBC7', '#9ACD32'))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]
  label_data <- go_enrich_df['category']
  color_assign <- c('BP' = Col1, 'CC' = Col2, 'MF' = Col3)
  circos.genomicTrack(
    plot_data, ylim = c(0, max(plot_data$rich.factor)), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')
      circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.05, 0.05), col = 'gray', lwd = 0.3) 
    } )
  circos.clear()
  category_legend <- Legend(
    labels = c("BP:biological process", "CC:cellular component","MF:molecular function" ),
    type = 'points', pch = NA, background = c(Col1,Col2,Col3), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background = c('#FDDBC7', '#9ACD32'), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
    title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(pvalue)')
  
  lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  pushViewport(viewport(x = 0.85, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  lgd_list_vertical <- packLegend(category_legend)
  pushViewport(viewport(x = 0.5, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  dev.copy(which = a) 
  dev.off()
  dev.off()
}

### go dotplot  ----------------------------------------------------------------
go_result <- go@result
# write.csv(go_result,file = "01.GO_result.csv",row.names = T)

### P值排序
BP <- go[which(go$ONTOLOGY=='BP'),]
CC <- go[which(go$ONTOLOGY=='CC'),]
MF <- go[which(go$ONTOLOGY=='MF'),]

### 将top5后的名字改为空值
go_BP  <- as.data.frame(BP)
go_BP [c(6:nrow(go_BP )),]$Description<-''
go_CC <- as.data.frame(CC)
go_CC [c(6:nrow(go_CC )),]$Description<-''
go_MF <- as.data.frame(MF)
go_MF [c(6:nrow(go_MF )),]$Description<-''

go_plot <- rbind(go_BP, go_CC, go_MF)
# ego_plot <- go_CC
# ego_plot <- go_BP
# ego_plot <- go_MF
colnames(go_plot)
go_plot$'-log10(p.adjust)' <- -log10(go_plot$p.adjust)
go_plot<-na.omit(go_plot)

if(T){
  p <- ggplot(go_plot,aes(x=Count,y=`-log10(p.adjust)`))+
    geom_point(aes(size=Count,color=ONTOLOGY),alpha=0.6)+
    scale_size(range=c(0.5,8))+
    scale_color_manual(values = c("#145390", "#68B69F", "#F08E59"))+ 
    # # scale_fill_gradient(low = "DeepSkyBlue",high = "Tomato")+
    # ylim(2,15)+
    # xlim(0,10)+
    theme_classic()+                              ###  也可以使用theme_bw()
    xlab("Count")+ylab("-log10(p.adjust)")+
    theme(axis.title.x =element_text(size=20, face = "bold"),
          axis.text.x =element_text(size=16, colour = "black"),
          axis.title.y =element_text(size=20, face = "bold"),
          axis.text.y=element_text(size=16, colour = "black"),
          strip.text = element_text(size = 16),
          plot.title=element_text(size=22),
          legend.title=element_text(size=16,face = "bold"),
          legend.text=element_text(size=14),
          
    )+
    theme(legend.position = c("right"))+
    geom_text_repel(fontface = "italic",
                    size =4,
                    color = "black",
                    data =go_plot,
                    aes(label = Description),
                    show.legend = FALSE )+
    facet_grid(.~ONTOLOGY)                        ### 图片分层
  p
  ggsave('03.GO_dotplot.pdf',p,w=16,h=10,bg='white')
  ggsave('03.GO_dotplot.png',p,w=16,h=10,units = "in", dpi = 600)
}


### go net ---------------------------------------------------------------------
### 双疾病，使用主要目标疾病
rt <- read.csv("../01_DEGs/02.DEGs_sig(GSE15605).csv", check.names=F, header=T, row.names = 1)
rt <- rt[gene.inter$symbol,]
rt$ID <- rownames(rt)
#colnames(rt)[1]<-"ID"
# colnames(rt)[2]<-"logFC"
rt<-rt[,c('ID','logFC')]

go.re <- ego@result
go.re <- mutate(go.re, 
                RichFactor = Count / as.numeric(sub("/\\d+","", BgRatio)))

go.re <- subset(go.re, pvalue < 0.05)
BP <- go.re[go.re$ONTOLOGY=='BP', ]
BP <- BP[order(BP$Count, decreasing = T),]
CC <- go.re[go.re$ONTOLOGY=='CC', ]
CC <- CC[order(CC$Count, decreasing = T),]
MF <- go.re[go.re$ONTOLOGY=='MF', ]
MF <- MF[order(MF$Count, decreasing = T),]


display_number = c(2, 2, 2)
ego_result_BP <- as.data.frame(BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(MF)[1:display_number[3], ]

go.re <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
colnames(go.re)
deg <- rt
colnames(deg) <- c('gene','logFC')

library(tidyverse)
temp1 <- go.re
temp1 <- dplyr::select(go.re, ID, geneID)
temp1 <- separate(temp1, geneID, into = c('var1','var2'), sep = '/')
ncol(temp1)
i <- 1
x <- 2
finally_data <- data.frame()
for (i in c(1 : nrow(temp1))) {
  for (x in c(2 : 3)) {
    temp2 <- temp1[i, ]
    temp3 <- data.frame(GO_ID = temp2[1, 1],
                        symbol = temp2[1, x])
    finally_data <- rbind(finally_data, temp3)
  }
  
}
finally_data <- na.omit(finally_data)

go.df <- merge(finally_data,go.re,by.x='GO_ID',by.y='ID')
go.df2 <- merge(deg,go.df,by.x='gene',by.y='symbol')

colnames(go.df2)
go.df2$abs_logFC <- abs(go.df2$logFC)
#treemap(go.df2,index=c("Description","gene"),vSize = "abs_logFC")
se_index <- c("Description","gene")
#root="GO Enrichment Terms"
library(tidyverse)
library(tidygraph)
library(ggraph)
# devtools::load_all('/data/nas2/Software/ccgraph/')
library(ccgraph)
library(treemap)
#devtools::install_github( "gaospecial/ccgraph")
nodes_GO <- gather_graph_node(go.df2,index = se_index, value = "abs_logFC",root="GO")
edges_GO <- gather_graph_edge(go.df2,index = se_index,root="GO")

head(nodes_GO,10)
head(edges_GO,10)

graph <- tbl_graph(nodes_GO,edges_GO)
library(ggsci)
if(T){
  gc <- ggraph(graph, layout = 'dendrogram', circular = TRUE) + 
    geom_edge_diagonal(aes(color = node1.node.branch,
                           filter=node1.node.level!="all"), 
                       alpha = 2/3,edge_width=0.5) + 
    geom_node_point(aes(size = node.size, 
                        alpha=0.5,
                        color = node.branch,
                        filter=node.level!="all")
    ) + 
    scale_size(range = c(0.5,20)) + 
    theme(legend.position = "none") + 
    
    scale_edge_color_brewer(palette = "Set1") + #
    scale_color_brewer(palette = "Set2") +
    
    geom_node_text(
      aes(
        x = 1.0175 * x, 
        y = 1.0175 * y, 
        label = node.short_name,
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        filter = leaf,
        color = node.branch
      ),
      size = 6, hjust = 'outward') +
    
    geom_node_text(
      aes(label=node.short_name,
          filter = !leaf & (node.level != "all"),
          color = node.branch),
      fontface="bold",
      #alpha=0.5,
      size=5,
      family="sans"
    ) + 
    theme(panel.background = element_rect(fill = NA)) +
    coord_cartesian(xlim=c(-1.25,1.25),ylim = c(-1.25,1.25))
  
  ggsave(filename = "04.GO_netplot.pdf", height = 10, width = 10, gc)
  ggsave(filename = "04.GO_netplot.png", height = 10, width = 10, gc, dpi = 600)
}


# 03.KEGG富集-------------------------------------------------------------------
KEGG_database <- 'hsa'
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 organism = KEGG_database, 
                 keyType = "kegg",
                 pAdjustMethod = 'none',
                 pvalueCutoff = 0.05
                 # ,qvalueCutoff = 0.05
)
save(kk, file = 'KEGG.Rdata')
load("KEGG.Rdata")
enrichKK <- DOSE::setReadable(kk, OrgDb= database,keyType='ENTREZID')
KEGG_results <- data.frame(enrichKK)
paste0("共富集到",nrow(KEGG_results),"条KEGG通路")
#  "共富集到63条KEGG通路"
write.csv(KEGG_results,'00.KEGG_result.csv', row.names = F)

### kegg dot -------------------------------------------------------------------
kegg.circ = 1
kegg.dot = 1
# kegg.chord = 1

kegg_plot<- KEGG_results
kegg_plot [c(6:nrow(kegg_plot )),]$Description<-''

kegg_plot$`-log10(p.adjust)` <- -log10(kegg_plot$p.adjust)
kegg_plot<-na.omit(kegg_plot)
kegg_plot$ONTOLOGY<-'KEGG'

if(kegg.dot == 1){
  p2 <- ggplot(kegg_plot,aes(x=`-log10(p.adjust)`,y= Description))+
    geom_point(aes(size=qvalue,color=ONTOLOGY),alpha=0.6)+
    scale_color_manual(values = c("#E17887")
                       # ,guide = FALSE    ### 可以去掉kegg legend
    ) + #柱状图填充颜色
    scale_size(range=c(.5,8))+
    # # scale_fill_gradient(low = "DeepSkyBlue",high = "Tomato")+
    # ylim(2,15)+
    # xlim(0,10)+
    theme_classic()+ xlab("-log10(p.value)")+ylab("Description")+
    theme(axis.title.x =element_text(size=20),
          axis.text.x =element_text(size=16, colour = "black"),
          axis.title.y =element_text(size=20),
          # axis.text.y=element_text(size=16),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 16),
          plot.title=element_text(size=22),
          legend.title=element_text(size=16),
          legend.text=element_text(size=14),
          #legend.position = "bottom",
          
    )+
    theme(legend.position = c("right"))+
    geom_text_repel(fontface = "italic",
                    size = 4,
                    color = "black",
                    data =kegg_plot,
                    aes(label = Description),
                    show.legend = FALSE )
  p2
  ggsave('05.KEGG_dotplot.pdf',p2,w=8,h=8)
  ggsave('05.KEGG_dotplot.png',p2,w=8,h=8,units = "in", dpi = 600)
}

### kegg circ ------------------------------------------------------------------
load(file = "KEGG.Rdata")
KEGG <- mutate(kk, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
kk <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_result <- kk@result
if(median(kegg_result$richFactor) < 0.1){kegg_result$richFactor<- kegg_result$richFactor*10}
hh <- as.data.frame(kegg_result)
rownames(hh) <- 1:nrow(hh)
hh <- hh[hh$pvalue <= 0.05,]
hh <- hh[order(hh$Count,decreasing = T),]

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
kk <- hh[c(1:5),] 

paste0(kk$Description,collapse = "（）；")
# "Asthma（）；Allograft rejection（）；Type I diabetes mellitus（）；Graft-versus-host disease（）；Intestinal immune network for IgA production"

if(kegg.circ == 1){
  kegg <-  data.frame(Category = "KEGG",ID = kk[,'ID'],Term = kk[,'Description'], 
                      Genes = gsub("/", ", ", kk[,'geneID']), adj_pval = kk[,'pvalue'])
  genelist <- data.frame(ID = rownames(select_DERNA), logFC = select_DERNA$`log2FoldChange`)
  row.names(genelist) = genelist[,1]
  library(GOplot)
  circ <- circle_dat(kegg, genelist)
  kk$pvalue <- -log10(kk$pvalue)
  go_enrich_df <- data.frame(
    category="KEGG", 
    gene_num.min = 0,
    gene_num.max = max(kegg_result$Count),
    gene_num.rich=kk$Count,
    "-log10.p" = kk$pvalue,
    up.regulated = 0,
    down.regulated = 0,
    rich.factor = kk$richFactor,
    Description=kk$ID)
  circ$sig <- ifelse(circ$logFC > 0,"up","down")
  tem1 <- go_enrich_df$Description
  for(i in 1:5){
    x <- tem1[i]
    circ1 <- circ[grep(x,circ$ID),]
    num1 <- nrow(circ1)
    tem2 <- as.data.frame(table(circ1$sig))
    if(nrow(tem2) == 2){
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[2]]
    }else if(tem2[1,1] == "up"){
      go_enrich_df[match(x,go_enrich_df$Description),6] = table(circ1$sig)[[1]]
    }else{
      go_enrich_df[match(x,go_enrich_df$Description),7] = table(circ1$sig)[[1]]}
  }
  
  go_enrich_df$Description <- factor(go_enrich_df$Description, levels = go_enrich_df$Description)
  rownames(go_enrich_df) <- go_enrich_df$Description
  library(ggplot2)
  pdf(file = paste0("06.KEGG_cicplot.pdf"),width = 13,height = 6)
  a <- dev.cur()   
  png(file = paste0("06.KEGG_cicplot.png"),width= 13, height= 6, units="in", res=300)
  dev.control("enable")
  par(family="Times")
  circle_size = unit(1, 'snpc')
  library(circlize)
  kegg.color <- sample(col1,1)
  circos.par(gap.degree = 2, start.degree = 90)
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max')]  #选择作图数据集，定义了 ko 区块的基因总数量范围
  ko_color <- c(rep(kegg.color, 10))  #定义分组颜色
  circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)  #一个总布局
  circos.track(
    ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  #圈图的高度、颜色等设置
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')  #ylim、xlim 用于指定 ko id 文字标签添加的合适坐标
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称
      circos.axis(h = 'top', labels.cex = 0.6, major.tick.percentage = 0.4, labels.niceFacing = FALSE)  #绘制外周的刻度线
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #将 ko id 文字标签添加在图中指定位置处
    } )
  colnames(go_enrich_df)[5] <- c("-log10(pvalue)")
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.rich', '-log10(pvalue)')]  #选择作图数据集，包括富集基因数量以及 p 值等信息
  label_data <- go_enrich_df['gene_num.rich']  #标签数据集，仅便于作图时添加相应的文字标识用
  p_max <- round(max(go_enrich_df$`-log10(pvalue)`)) + 1  #定义一个 p 值的极值，以方便后续作图
  colorsChoice <- colorRampPalette(c('white','#4393C3'))  #这两句用于定义 p 值的渐变颜色
  color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))
  
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #区块的长度反映了富集基因的数量，颜色与 p 值有关
      ylim = get.cell.meta.data('ycenter')  #同上文，ylim、xlim、sector.name 等用于指定文字标签（富集基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),1]
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #将文字标签添（富集基因数量）加在图中指定位置处
    } )
  go_enrich_df$all.regulated <- go_enrich_df$up.regulated + go_enrich_df$down.regulated
  go_enrich_df$up.proportion <- go_enrich_df$up.regulated / go_enrich_df$all.regulated
  go_enrich_df$down.proportion <- go_enrich_df$down.regulated / go_enrich_df$all.regulated
  
  go_enrich_df$up <- go_enrich_df$up.proportion * go_enrich_df$gene_num.max
  plot_data_up <- go_enrich_df[c('Description', 'gene_num.min', 'up')]
  names(plot_data_up) <- c('Description', 'start', 'end')
  plot_data_up$type <- 1  #分配 1 指代上调基因
  
  go_enrich_df$down <- go_enrich_df$down.proportion * go_enrich_df$gene_num.max + go_enrich_df$up
  plot_data_down <- go_enrich_df[c('Description', 'up', 'down')]
  names(plot_data_down) <- c('Description', 'start', 'end')
  plot_data_down$type <- 2  #分配 2 指代下调基因
  
  plot_data <- rbind(plot_data_up, plot_data_down)
  label_data <- go_enrich_df[c('up', 'down', 'up.regulated', 'down.regulated')]
  color_assign <- colorRamp2(breaks = c(1, 2), col = c('#FDDBC7', '#9ACD32'))
  circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  #这里紫色代表上调基因，蓝色代表下调基因，区块的长度反映了上下调基因的相对占比
      ylim = get.cell.meta.data('cell.bottom.radius') - 0.5  #同上文，ylim、xlim、sector.name 等用于指定文字标签（上调基因数量）添加的合适坐标
      xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),3]
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #将文字标签（上调基因数量）添加在图中指定位置处
      xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
      sector.name = label_data[get.cell.meta.data('sector.index'),4]
      circos.text(xlim, ylim, sector.name, cex = 0.6, niceFacing = FALSE)  #类似的操作，将下调基因数量的标签也添加在图中
    } )
  plot_data <- go_enrich_df[c('Description', 'gene_num.min', 'gene_num.max', 'rich.factor')]  #选择作图数据集，标准化后的富集得分
  label_data <- go_enrich_df['category']  #将通路的分类信息提取出，和下一句一起，便于作图时按分组分配颜色
  color_assign <- c('KEGG' = kegg.color)
  circos.genomicTrack(
    plot_data, ylim = c(0, max(plot_data$rich.factor)), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  #圈图的高度、颜色等设置
    panel.fun = function(region, value, ...) {
      sector.name = get.cell.meta.data('sector.index')  #sector.name 用于提取 ko id 名称，并添加在下一句中匹配 ko 对应的高级分类，以分配颜色
      circos.genomicRect(region, value,col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  #绘制矩形区块，高度代表富集得分，颜色代表 ko 的分类
      circos.lines(c(0, max(region)), c(0.05, 0.05), col = 'gray', lwd = 0.3)  #可选在富集得分等于 0.5 的位置处添加一个灰线
    })
  circos.clear()
  category_legend <- Legend(
    labels = c("KEGG Pathway"),
    type = 'points', pch = NA, background = c(kegg.color), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  updown_legend <- Legend(
    labels = c('Up-regulated', 'Down-regulated'), 
    type = 'points', pch = NA, background =c('#FDDBC7', '#9ACD32'), 
    labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
  pvalue_legend <- Legend(
    col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                         colorRampPalette(c('white','#4393C3'))(6)),
    legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
    title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(pvalue)')
  
  lgd_list_vertical <- packLegend(updown_legend, pvalue_legend)
  pushViewport(viewport(x = 0.85, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  lgd_list_vertical <- packLegend(category_legend)
  pushViewport(viewport(x = 0.5, y = 0.5))
  grid.draw(lgd_list_vertical)
  upViewport()
  dev.copy(which = a) 
  dev.off()
  dev.off()
}
