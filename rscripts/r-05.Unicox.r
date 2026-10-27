rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "05_Unicox"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)

pacman::p_load(survival, survminer, tidyr, dplyr, forestplot)
# 01 获取数据集  ---------------------------------------------------------------
## 表达矩阵
dat <- read.csv('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/SKCM/01.dat.fpkm.csv', row.names = 1)
# dat <- log2(dat+.5)
colnames(dat) <- gsub(".","-",colnames(dat),fixed = T)

## 生存信息
survival <- read.csv("../00_rawdata/SKCM/04.survival.csv")
dat <- dat[,which(colnames(dat)%in%survival$sample)]

##
gene <- read.csv("../04_MR/02.MR_gene.csv")
colnames(gene)

train_dat <- dat
train_dat <- t(train_dat)

## 合并生存数据
train_dat <- t(dat[gene$SYMBOL, colnames(dat) %in% survival$sample]) %>% data.frame()
train_dat$sample <- rownames(train_dat)
train_dat <- merge(survival,train_dat,by='sample')
rownames(train_dat) <- train_dat$sample
train_dat <- train_dat[,c(-1)]

colnames(train_dat)

# 02.单因素cox -----------------------------------------------------------------
colnames_sum <- colnames(train_dat)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
# colnames_sum <- gsub("RGS5.1","RGS5_1",colnames_sum)
colnames(train_dat) <- colnames_sum
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]

# Surv()函数产生一个生存对象  生存时间对生存的影响 对每一个变量构建生存分析公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))  #as.formula(). 将字符串转换成公式。构建formula对象
# coxph函数用于计算cox模型 循环对每一个特征做cox回归分析
univ_models <- lapply(univ_formulas,
                      function(x) {coxph(x, data = train_dat)})


## 02.1 COX建模并进行PH假定检验 ------------------------------------------------
univ_cox_zph <- lapply(univ_models,
                       function(x) {
                         y <- cox.zph(x)
                         #获取p值
                         p<-y$table[-nrow(y$table),][3]
                         chisq <- y$table[-nrow(y$table),][1]
                         df <- y$table[-nrow(y$table),][2]
                         res<-c(chisq,df,p)
                         names(res)<-c("chisq","df","p")
                         return(res)
                       })
univ_cox_zph <- t(as.data.frame(univ_cox_zph, check.names = FALSE))
univ_cox_zph <- as.data.frame(univ_cox_zph)
univ_cox_zph <- rownames(univ_cox_zph)[univ_cox_zph$p>0.05]

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <- signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

## coef是公式中的回归系数b（有时候也叫beta值）。exp(coef)是cox模型中的风险比（HR）
## z代表wald统计量，是coef除以其标准误se(coef)。ower .95 upper .95则是exp(coef)的95%置信区间，可信区间越窄，可信度越高，你的实验越精确，越是真理。
res_mod <- t(as.data.frame(univ_results, check.names = FALSE))
res_mod <- as.data.frame(res_mod)
res_mod <- res_mod[rownames(res_mod)%in%univ_cox_zph, ]
res_results_0.05 <- res_mod[which(as.numeric(res_mod$p.value) < 0.05),]
res_results_0.05 <- res_results_0.05[-6,]
res_results_0.05 <- na.omit(res_results_0.05)
dim(res_results_0.05)
paste0("共", nrow(res_results_0.05),"个基因通过单因素分析。") 
# 共2个基因通过单因素分析。

write.csv(res_results_0.05, file = "01.univariate_cox_result_0.05.csv")


##  选取0.05结果,HR不等于1的进行后续分析
library(tidyr)
res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))

##  选取HR前15基因进行图片绘制
# res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$p.value),]
# res_results_0.05_2 <- res_results_0.05_2[1:15,]
hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
# view the last line
nrow(tabletext) + 1




fores.plot <- forestplot(labeltext=tabletext, 
                         graph.pos=4,  #为Pvalue箱线图所在的位置
                         is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
                         col=fpColors(box="#ed2f6a", lines="#2177b8", zero = "gray"),
                         title = "Cox ProportionalHazards Model",
                         mean=c(NA,NA,res_results_0.05_2$HR),
                         lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
                         upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
                         boxsize=0.15,lwd.ci=2,   #箱子大小，线的宽度
                         ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
                         zero=1,lwd.zero=2,    #zero线宽 基准线的位置
                         colgap=unit(5,"mm"),    #列间隙
                         xticks = c(0,0.5,1.0,1.5,2.0),
                         # xticks = c(0, 1,2,5,10,20,30), #横坐标刻度
                         lwd.xaxis=2,            #X轴线宽
                         lineheight = unit(0.8,"cm"), #固定行高
                         graphwidth = unit(.5,"npc"), #图在表中的宽度比例
                         cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                         hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2),
                                           "9" = gpar(col = "black", lty = 1, lwd = 2)),
                         # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                         #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                         #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
                         # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
                         #fpTxtGp函数中的cex参数设置各个组件的大小
                         txt_gp=fpTxtGp(label=gpar(cex=1),
                                        ticks=gpar(cex=0.8, fontface = "bold"),
                                        xlab=gpar(cex = 1, fontface = "bold"),
                                        title=gpar(cex = 1, fontface = "bold")),
                         xlab="Hazard Ratio",
                         grid = T) # 垂直于x轴的网格线，对应每个刻度

pdf(file = "02.univariate_cox_forest.pdf", height = 6, width = 10, onefile = F)
print(fores.plot)
dev.off()

png(filename = "02.univariate_cox_forest.png", height = 6, width = 10, units = "in", res = 600)
print(fores.plot)
dev.off()

# 提取通过单因素的基因 ----------------
hubgene <- data.frame(symbol = rownames(res_results_0.05))
write.csv(hubgene, "03.hubgene.csv", row.names = F)
