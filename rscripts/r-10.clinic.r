rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "10_clinical"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)

###---- load packages ----####
library(tidyr); library(ggplot2); library(tidyverse); library(ggsci); library(ggpubr)
library(gghalves); library(ggthemes); library(patchwork); library(lance); library(cowplot)
library(survminer); library(survival)

###---- data set ----####
work.path <- "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/10_clinical/"
single.plot <- file.path(work.path, "singleplot")
km.single.plot <- file.path(work.path, "km_singleplot")

if( !dir.exists(single.plot)){dir.create(single.plot)}
if( !dir.exists(km.single.plot)){dir.create(km.single.plot)}

###---- load clinic data ----####
clinic.dat <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/09_independent/01.phenotype.csv")
str(clinic.dat)

clinic.dat <- lapply(clinic.dat, function(x) if(is.character(x)) as.factor(x) else(x))
clinic.clear <- clinic.dat %>% as.data.frame()
str(clinic.clear)
clinic.clear$risk <- factor(clinic.clear$risk,levels = c('Low', "High"))

###---- wilcox test of risk group ----####
train_phenotype <- clinic.clear
colnames(train_phenotype)


# create wilcox group
create_paired_group <- function(members) {
  # 根据所有组成员创建分组list
  # param1: members 所有成员向量，类似于c("A", "B", "C")
  # return: paired_group_list 配对成员分组的list，类似于list(c("A", "B"), c("A", "C"), c("B", "C"))
  paired_group_list <- list()
  end <- length(members)
  for (i in 1:(end - 1)) {
    for (j in (i + 1):end) {
      paired_members <- c(members[i], members[j])
      paired_group_list[[length(paired_group_list) + 1]] <- paired_members
    }
  }
  
  return(paired_group_list)
}
# 提取变量名称
variable <- function(data){
  var <- unique(na.omit(data))
  var <- var[order(var[,1]),] %>% as.character()
  return(var)
}

# wilcox test plot
wil.plot <- function(wil.dat, my_comparisons, title ,mytheme){
  if (length(my_comparisons) <= 1) {
    p <- ggplot(wil.dat, aes(x = group, y = riskScore, fill = group)) +
      geom_half_violin(trim = F, color = NA, side = "r") +  
      #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
      #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
      geom_half_boxplot(width=0.4,position=position_dodge(0.9),outlier.shape = NA,fill='white',side = "r",
                        errorbar.draw = F)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
      geom_half_point_panel(side = "l",shape = 21,color = "white",size = 3)+ 
      scale_fill_manual(values = c("#CED0E7","#D5B8AC","#B0ACAD","#7D7667","#647167"), name = "Group")+
      labs(title="", x="", y = "", size = 20) +
      geom_signif(comparisons = my_comparisons,test = wilcox.test,map_signif_level = T,
                  y_position = max(wil.dat$riskScore) +.6)+
      ggtitle(paste0(title))+
      theme_bw()+ mytheme 
    # +guides(fill='none')
  }else{
    p <- ggplot(wil.dat, aes(x = group, y = riskScore, fill = group)) +
      geom_half_violin(trim = F,color = NA,side = "r") + 
      #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
      #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
      geom_half_boxplot(width = 0.4, position=position_dodge(0.9), outlier.shape = NA, fill='white', side = "r",
                        errorbar.draw = F)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
      geom_half_point_panel(side = "l", shape = 21, color = "white", size = 3)+
      scale_fill_manual(values= c("#CED0E7","#D5B8AC","#B0ACAD","#7D7667","#647167"), name = "Group")+
      labs(title="", x="", y = "", size=20) +
      geom_signif(comparisons = my_comparisons, test = wilcox.test, map_signif_level = T, step_increase = 0.1)+
      ggtitle(paste0(title))+
      theme_bw()+ mytheme
    # +guides(fill='none')
  }
  return(p)
}

# mytheme
mytheme <- theme(plot.title = element_text(hjust = 0.5, colour = "black",face = "bold", size = 18),
                 axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", face = "bold", size = 12),
                 axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size=12),
                 axis.title.x = element_text(size = 16, face = "bold"),
                 axis.title.y = element_text(size = 16, face = "bold"),
                 legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 12),
                 legend.title = element_text(face = "bold", size = 12),
                 legend.position = "top",
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())

# extract plot clinic
colnames(train_phenotype)
clinical.name <- colnames(train_phenotype)[which(!colnames(train_phenotype) %in% c("id", "riskScore", "risk", "OS", "OS.time", "gender", "group", "age"))]

# plot list 
plot.list <- list()
# i <- "tumor.stage"
j <- 1
for (i in clinical.name) {
  my_comparisons <- create_paired_group(variable(train_phenotype[i]))
  level <- variable(train_phenotype[i])
  plot.dat <- data.frame(riskScore = train_phenotype$riskScore,
                         group = factor(train_phenotype[,i], levels = level))
  plot.dat <- na.omit(plot.dat)
  p <- wil.plot(plot.dat, my_comparisons, i, mytheme)
  ggsave(file.path(single.plot, paste0("0", j, ".", i, ".pdf")), width = 8, height = 8)
  ggsave(file.path(single.plot, paste0("0", j, ".", i, ".png")), width = 8, height = 8, dpi = 600)
  plot.list[[i]] <- p
  j <- j+1
}

combined.plot <- plot_grid(plotlist = plot.list, ncol = 2)
combined.plot

ggsave(filename = '01.clinical_box.pdf', combined.plot, width = 14, height = 14)
ggsave(filename = '01.clinical_box.png', combined.plot, width = 14, height = 14, dpi = 600, bg = "white")


###---- kaplan meier of clinic ----####
km.clil.data <- clinic.clear

km.clil.data$risk <- factor(km.clil.data$risk, levels = c('High','Low'))
clinical.name


j <- 1
km.plot.list <- list()
for (i in clinical.name) {
  var.name <- variable(km.clil.data[clinical.name[j]])
  tmp.clil.name <- clinical.name[j]
  for (n in c(1:length(var.name))) {
    tmp.var <- var.name[n]
    # if(tmp.var == "tumor.stage4"){break}
    # if(tmp.var == "histologic.stage1"){next}
    tmp.dat <- km.clil.data[which(km.clil.data[tmp.clil.name] == tmp.var), ]
    tmp.kmfit <- survfit(Surv(OS.time,OS) ~ risk, data = tmp.dat)
    p <- ggsurvplot(tmp.kmfit,
                    pval = TRUE, 
                    conf.int = F,
                    legend.labs = c("High risk","Low risk"),
                    legend.title = "risk group",
                    title = paste0("Kaplan Meier of ",tmp.clil.name," ",tmp.var),
                    font.main = c(15,"bold"),
                    risk.table = TRUE, 
                    risk.table.col = "strata",  
                    linetype = "strata", 
                    surv.median.line = "hv", 
                    ggtheme = theme_bw(), 
                    palette = c("#A73030FF", "#0073C2FF")) 
    ggsave(file.path(km.single.plot, paste0(tmp.clil.name,"_", n, ".pdf")), p$plot, width = 8, height = 8)
    ggsave(file.path(km.single.plot, paste0(tmp.clil.name,"_", n, ".png")), p$plot,width = 8, height = 8, dpi = 600)
    km.plot.list[[tmp.var]] <- p
  }
  j <- j + 1
}

km.plot.list.plot <- list()
km.name <- names(km.plot.list)
for (i in c(1:length(km.name))) {
  km.plot.list.plot[[i]] <- km.plot.list[[i]]$plot
}

combined.plot.km <- plot_grid(plotlist = km.plot.list.plot, ncol = 3)
combined.plot.km

ggsave(filename = '02.clinical_km.pdf', combined.plot.km, width = 22, height = 16)
ggsave(filename = '02.clinical_km.png', combined.plot.km, width = 22, height = 16, dpi = 600, bg = "white")

