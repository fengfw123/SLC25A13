rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "12_cibersort"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)

###---- load packages ----####
library(ggpubr); library(IOBR); library(magrittr); library(RColorBrewer); library(rstatix); library(reshape2)
library(corrplot); library(psych); library(Hmisc); library(ggplot2); library(ggcor); library(ggpubr)

###---- load data ----####
train.data <- read.csv(file = "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/SKCM/01.dat.fpkm.csv", header=TRUE, row.names=1, check.names = FALSE)
group <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/07_Riskmodel/01.train_riskScore.csv", row.names = 1) %>% dplyr::select(c("id", "risk"))
colnames(group) <- c('id', 'group')

dat <- train.data[,group$id]
exprSet <- dat

cibersort <- deconvo_tme(eset = exprSet, method = "cibersort", arrays = T, perm = 200 )

cibersort.res <- cibersort %>% as.data.frame()
colnames(cibersort.res) %<>% gsub("_CIBERSORT","",.)
colnames(cibersort.res) %<>% gsub("_"," ",.)
rownames(cibersort.res) <- cibersort.res$ID
cibersort.res <- cibersort.res[,-c(1,24:26)]

### 输出结果文件
write.csv(cibersort.res,'01.CIBERSORT_result.csv')
save.image("result.RData")
load("result.RData")


###---- cibersort stacked diagrams ----####
stack.dat <- t(cibersort.res) %>% as.data.frame()
stack.dat$cell_type <- rownames(stack.dat)

mypalette <- colorRampPalette(brewer.pal(8,"Accent"))

# cibersort stacked
stack.dat %>% gather(id, fraction, -cell_type) %>%  merge(group, by = 'id') %>%
  ggplot(aes(x = id, y = fraction, fill = cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x = '', y = 'Relative Percent', fill = '')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(~ group, scales = "free", space = "free")

# extract plot
save.p <- recordPlot()
pdf('02.cibersort_stack.pdf',width = 10,height = 8, family = 'Times')
replayPlot(save.p)
dev.off()

png('02.cibersort_stack.png',w=10,h=8,units='in',res=600, family = 'Times')
replayPlot(save.p)
dev.off()

###---- cibersort wilcox test ----####
# delete < 0
wil.dat <- cibersort.res[,which(colSums(cibersort.res) > 0)]
# order
wil.dat <- wil.dat[group$id,order(colnames(wil.dat),decreasing = F)]
identical(group$id, rownames(wil.dat))

# merge data
wil.dat <- wil.dat %>% as.data.frame %>% tibble::rownames_to_column(var = "id")
identical(wil.dat$id, group$id)

# 不做基因相关性，需去掉基因
wil.dat <- merge(wil.dat, group, by = "id")

# wilcox test anlysis
dat.test <- tidyr::gather(wil.dat, Cell, `Cell composition`, -c(id, group))
colnames(dat.test)
stat.res <- dat.test %>% 
  group_by(Cell) %>% 
  wilcox_test(`Cell composition` ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat.res  

# ouput result
write.csv(stat.res, file = "03.wilcoxon.csv", row.names = F, quote = F)

# extract p < 0.05 cells
dif.cell <- stat_res[stat.res$p<0.05, ]
paste0("共获得 ", nrow(dif.cell)," 个差异免疫细胞")
# "共获得 8 个差异免疫细胞"


###---- cibersort wilcox test box plot ----####
# organize data
box.dat <- dat.test[dat.test$Cell %in% dif.cell$Cell, ]
box.cell <- dif.cell

# box plot
if(T){
  p <- ggboxplot(box.dat, x = "Cell", y = "Cell composition", fill = "group",
                 palette = c("#69bcce", "#fed53b"), outlier.size = 0.8, bxp.errorbar = T) +
    ylim(0,.6)+
    stat_pvalue_manual(box.cell , x = "Cell", y.position = .6, size = 3.5, color = "black",
                       family = "Times", label = "p.signif") +
    theme_bw()+
    labs(title = "Differential immune cells", x = "", y = "Cell composition", color = "") +
    theme(plot.title = element_text(hjust =0.5, colour = "black", face = "bold", size = 16),
          axis.text.x = element_text(angle = 0, hjust = 0.5, colour ="black", face = "bold", size = 12), 
          axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12), 
          axis.title.x = element_text(size = 16,face = "bold"),
          axis.title.y = element_text(size = 16,face = "bold"),
          legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 14),
          legend.title = element_text(face = "bold", size = 14),
          text = element_text(family = 'Times'),
          legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p
}

# ouput box plot
ggsave(filename = "04.differential_immune_cells.pdf", p, width = 14, height = 7)
ggsave(filename = "04.differential_immune_cells.png", p, width = 14, height = 7, dpi = 600, bg = "white")


###---- cibersort cor for riskscore and diffcell ----####
# riskscore
riskscore <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/07_Riskmodel/01.train_riskScore.csv", row.names = 1) %>% dplyr::select(c("id", "riskScore"))
identical(riskscore$id, rownames(cor.dat))
# organize cor data
cor.dat <- wil.dat %>% column_to_rownames(var = 'id')

# extract diffcell
cor.dat <- cor.dat[ ,dif.cell$Cell]
cor.dat$id <- rownames(cor.dat)

# merge cell data and riskscore
cor.dat <- merge(cor.dat, riskscore, by = "id", all = TRUE)
rownames(cor.dat) <- cor.dat$id
cor.dat <- cor.dat[, colnames(cor.dat) != "id"]

#  cor of spearman or pearson
cor.res <- cor(cor.dat, method = "spearman")
corp <- cor_pmat(cor.dat)
write.csv(cor.res, '05.riskcor.csv', quote = F)

if(T){
  geng_num <- 1
  cell_num <- nrow(dif.cell)
  mycolors <- c( "#104680", "#59829e", "white", "#f1939c", "#d11a2d") 
  corr_plot <- quickcor(cor.res,
                        # method = "spearman",                     # 默认方法为 pearson
                        is.corr = T,                           # 输入得数据是否是相关性系数结果，默认F，表示输入得是元素数据；
                        type = c("lower", "upper", "full")[1], # 修改中括号内数值1~3选择不同类型的相关性矩阵；
                        # p.mat = corp,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(100),
                        tl.col = 'black',
                        tl.offset = 0.4,
                        number.font = 2,
                        cl.cex = 0.8,
                        number.cex = 1.2
  )+
    # 可在下面任选一个图形，其他的在代码最前面加#号
    geom_square(colour = "#5a5c5f", size = 0.6, linetype = 2) +  #绘制正方形热图，描边粗细、线条类
    # geom_circle2(colour = "white",size = 0.25, linetype = 1)+   #绘制圆形热图 ，描边粗细、线条类型
    # geom_ellipse2(colour = "#5a5c5f", size = 0.6, linetype = 2)+ #椭圆形热图，描边粗细、线条类型
    # geom_colour(size = 0.3, linetype = 1)+                       #经典热图，描边粗细、线条类型
    # geom_pie2()+ #绘制饼图
    # geom_star(colour = "white",size = 0.25)+ #绘制五角星热图
    geom_number(aes(num = r), data = get_data(show.diag = FALSE, abs(r) > 0.3), size = 3.5, color = "white")+ # 添加r值文字
    scale_fill_gradientn(colours = mycolors, limits = c(-1,1))+  #自定义颜色
    # breaks = c(-1, 0, 1) +  #自定义刻度值
    labs(title = "Gene and Cell Spearman correlation",
         fill = "R(Spearman)")+
    theme(
      plot.title = element_text(color = "black", size = 15 ,face = "bold", hjust = .5),
      axis.text.x = element_text(color = c(rep(mycolors[1], cell_num), rep(mycolors[length(mycolors)], geng_num)), # X轴标签颜色
                                 #face = "bold",
                                 size = 12),                       # 设置加粗字体和大小
      axis.text.y = element_text(color = c(rep(mycolors[length(mycolors)], geng_num), rep(mycolors[1], cell_num)), # Y轴标签颜色
                                 #face = "bold",
                                 size = 12),                       # 设置加粗字体和大小
      legend.title = element_text(size = 13, color = "black", face = "bold"),   # 图例标题
      legend.text = element_text(size = 12, color = "black"),    
      legend.position = "right"                                    # 图例位置
    )
  corr_plot
  ggsave("06.risk_corheatmap.pdf", corr_plot, height = 14, width = 14)
  ggsave("06.risk_corheatmap.png", corr_plot, height = 14, width = 14, dpi = 600, units = "in")
}


