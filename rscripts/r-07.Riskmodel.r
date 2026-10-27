rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "07_Riskmodel"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)

###---- load package ----####
library(ggplot2); library(ggthemes); library(Ipaper); library(survival); library(survminer)
library(dplyr); library(survivalROC); library(ComplexHeatmap)

###---- load rawdata ----#### 
# model riskscore
RS.train <- read.csv('../06_Model101/Results/02.RS_train.csv',row.names = 1)

# survival data
survival <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/SKCM/04.survival.csv")
rownames(survival) <- survival$sample

# risk data
survival <- survival[rownames(RS.train), ]
identical(rownames(RS.train), rownames(survival))
risk <- merge(survival, RS.train, by.x = "row.names", by.y = "row.names")
risk$riskScore <- risk$StepCox.backward..RSF
risk <- risk[,c("Row.names","OS","OS.time","riskScore")]
colnames(risk)[1] <- 'id'

risk <- risk[risk$OS.time > 60, ]
# cut point
# res.cut <- surv_cutpoint(risk_out2,time = 'OS.time',event = 'OS',variables = 'riskScore_out',minprop = 0.46)
# res.cut
# cutpoint <- res.cut$cutpoint[1]
# ###截断值计算
# risk_out2$risk_out <- ifelse(risk_out2$riskScore>cutpoint$cutpoint,'High Risk','Low Risk')

# median 
risk$risk <- ifelse(risk$riskScore > median(risk$riskScore), 'High', 'Low')
write.csv(risk, file = "01.train_riskScore.csv", quote = F)
table(risk$risk)


###---- train kaplan meier ----####
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data =  risk)
surv_pvalue(kmfit)

train_survival <- ggsurvplot(kmfit,
                            pval = TRUE, 
                            conf.int = F,
                            legend.labs=c("High risk","Low risk" ),
                            legend.title="Risk score",
                            title="Test Kaplan Meier",
                            font.main = c(15,"bold"),
                            risk.table = TRUE, 
                            risk.table.col = "strata", 
                            linetype = "strata", 
                            surv.median.line = "hv", 
                            ggtheme = theme_bw(), 
                            palette = c("#A73030FF", "#0073C2FF"))
train_survival$table <- train_survival$table +
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", color = 'black',family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
train_survival$plot <- train_survival$plot +
  labs(x = "Time (days)") +
  labs(y = "Test survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", color = 'black', family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(train_survival)

pdf('01.train_KaplanMeier.pdf', width = 8, height = 8, onefile = F, family = 'Times')
print(train_survival)
dev.off()
png('01.train_KaplanMeier.png', width = 8, height = 8, units = "in", family = 'Times', res = 600)
print(train_survival)
dev.off()

###---- ROC curve ----####
ROC <- risk

# roc function 
multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore,
                           predict.time=single_time,
                           method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'= rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC(time_vector = c(365*seq(3,7,2)), risk_score_table = risk)
table(for_multi_ROC$AUC)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# calculate auc fuction 
auc.value <- function(data, time){
  auc <- round(data[which(data$Time == time),5][1],2)
}

# data set 
time <- c(365*seq(3,7,2))
year <- c(3,5,7)
auc <- c(auc.value(for_multi_ROC, 365*3), auc.value(for_multi_ROC, 365*5), auc.value(for_multi_ROC, 365*7))
roc.list <- list()
color <- c("#BC3C29FF", "#0072B5FF", "#E18727FF")

# extract roc plot
for (i in 1:length(auc)) {
  roc.list[[i]] <- ggplot(for_multi_ROC[for_multi_ROC$Time == time[i],],
                          aes(x = False_positive,
                              y = True_positive,
                              label = Cut_values,
                              color = Time)) +
    scale_color_manual(breaks = time[i],
                       labels = paste0(year[i], " years"),
                       values = color[i]) +
    plotROC::geom_roc(labels = F, stat = 'identity') +
    style_roc() +
    geom_abline(slope = 1, intercept = 0, color = 'gray', linetype = 2) +
    theme_classic() +
    labs(title = "train ROC Method : Kaplan Meier", 
         subtitle = paste0("Years of ", year[i])) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size=12, color = "black", family = "Times", face = "bold"),
          axis.text.y = element_text(size=12, color = "black", family = "Times", face = "bold"),
          axis.title.x = element_text(size=14, color = "black", family = "Times", face  = "bold"),
          axis.title.y = element_text(size=14, color = "black", family = "Times", face = "bold"),
          legend.title = element_text(size=12, color = "black", family = "Times", face = "bold"),
          legend.text = element_text(size = 11, color = "black", family = "Times", face = "bold"),
          plot.title = element_text(size = 15, hjust = 0.5, color = "black", face = "bold"),
          plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black", face = "bold")) +
    annotate("text", x = 0.75, y = c(0.25),
             label = c(paste0('AUC of ', year[i], ' years = ', auc[i])))
  
}

roc.list[[1]]
ROC <- wrap_plots(roc.list, nrow = 1)

ggsave('02.train_ROC.pdf', ROC, width = 16, height = 6)
ggsave('02.train_ROC.png', ROC, width = 16, height = 6, dpi = 600)


###---- survival state plot ----####
rt <- risk
rt <- rt[order(rt$riskScore),]
riskClass <- rt[,"risk"]
lowLength <- length(riskClass[riskClass=="Low"])
highLength <- length(riskClass[riskClass=="High"])
color <- as.vector(rt$OS)
color <- ifelse(color == 1, "#dc143c", "#6495ed")
# color[color==1] <- "#dc143c"     ##1:dead, 0:alive
# color[color==0] <- "#6495ed"
table(color)

pdf(file = "03.train_survStat.pdf", width = 10, height = 6)
par(pin = c(4,4), mar = c(6,6,6,1),family = 'serif')
# pin: 以英寸表示的图形尺寸（宽和高）
# mai: 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英寸
# mar: 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英分。默认（5,4,4,2）+ 0.1
plot(rt$OS.time/365,
     pch = 19,
     xlab = "Patients (increasing risk socre)",
     ylab = "Time (years)",
     col = color,
     cex.axis = 1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab = 2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main = 2.0,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     #main="\n Risk Score Distribution"
)
abline(v = lowLength,lty = 2)
legend("topright", c(paste("Dead"),
                     paste("Alive")),
       x.intersp = 1, y.intersp = 0.8,
       pch = 20, col = c("#dc143c","#6495ed"),
       bty = "1",# bty框的类型  ,n:无边框  1：有边框
       seg.len = 1,cex = 1.5,
       text.font = 2)
title(main = "train survival state", cex.main = 2)
dev.off()


png(file = "03.train_survStat.png",width = 10,height = 6,units='in', res = 600)
par(pin = c(4,4), mar = c(6,6,6,1),family = 'serif')
plot(rt$OS.time/365,
     pch = 19,
     xlab = "Patients (increasing risk socre)",
     ylab = "Time (years)",
     col = color,
     cex.axis = 1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab = 2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main = 2.0,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     #main="\n Risk Score Distribution"
)
abline(v = lowLength,lty = 2)
legend("topright", c(paste("Dead"),
                     paste("Alive")),
       x.intersp = 1, y.intersp = 0.8,
       pch = 20, col = c("#dc143c","#6495ed"),
       bty = "1",# bty框的类型  ,n:无边框  1：有边框
       seg.len = 1,cex = 1.5,
       text.font = 2)
title(main = "train survival state", cex.main = 2)
dev.off()

###---- survival curve ----####
rt <- risk
rt <- rt[order(rt$riskScore),]
riskClass <- rt[,"risk"]
lowLength <- length(riskClass[riskClass == "Low"])
highLength <- length(riskClass[riskClass == "High"])
line <- rt[,"riskScore"]


pdf(file = "04.train_survcurve.pdf", width = 10, height = 6)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(line,
     type = "p",
     pch = 20,
     xlab = "Patients (increasing risk socre)",
     ylab = "Risk score",
     col = c(rep("#6495ed",lowLength),
           rep("#dc143c",highLength)),
     cex.axis = 1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab = 2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main = 2.0,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2)
legend("topleft",c(paste("High"),
                   paste("Low")),
       x.intersp = 1, y.intersp = 0.8,
       pch = 20,col=c("#dc143c", "#6495ed"),
       bty = 1,# bty框的类型
       seg.len = 1,cex = 1.5,
       text.font = 2)#
abline(h = median(rt$riskScore), v = lowLength, lty = 2)
title(main = "train survival curve", cex.main = 2)
dev.off()

png(file="04.train_survcurve.png",width = 10,height = 6,units = 'in', res = 600 )
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(line,
     type = "p",
     pch = 20,
     xlab = "Patients (increasing risk socre)",
     ylab = "Risk score",
     col = c(rep("#6495ed",lowLength),
             rep("#dc143c",highLength)),
     cex.axis = 1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab = 2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main = 2.0,   ##标题的缩放倍数。
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub = 2)
legend("topleft",c(paste("High"),
                   paste("Low")),
       x.intersp = 1, y.intersp = 0.8,
       pch = 20,col=c("#dc143c", "#6495ed"),
       bty = 1,# bty框的类型
       seg.len = 1,cex = 1.5,
       text.font = 2)#
abline(h = median(rt$riskScore), v = lowLength, lty = 2)
title(main = "train survival curve", cex.main = 2)
dev.off()



###---- gene expression ----####
dat <- read.csv(file = "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/SKCM/01.dat.fpkm.csv", header=TRUE, row.names=1, check.names = FALSE)
colnames(dat) <- gsub('.', '-', colnames(dat), fixed = T)

gene <- read.csv(file = "../06_Model101/Results/03.fea_df.csv")
gene <- gene[gene$algorithm == "StepCox[backward]+RSF",]$features


rdt <- risk
rdt <- rdt[order(risk$risk,decreasing = T),]
expr <- subset(dat, select = rdt$id)
a <- data.frame(row.names = risk$id, Risk = risk$risk)

select.expr <- expr[gene,]
ann.color <- list(Risk = c(High = "#dc143c",Low ="#6495ed"))
matrix <- t(scale(t(select.expr)))#归一化
matrix[matrix < -1] <- -1
matrix[matrix > 1] <- 1

pheatmap::pheatmap(matrix,main = 'Train data gene expression',
                   fontfamily = "Times",
                   fontsize = 10,
                   annotation_colors = ann.color,
                   annotation_col = a,
                   annotation_names_col = F,
                   annotation_names_row = F,
                   color = colorRampPalette(c("navy", "white","firebrick1"))(500),
                   cluster_rows = F,
                   cluster_cols = F,
                   show_colnames = F,
                   width = 6,height = 4,
                   file = '05.heatmap.pdf')
pheatmap::pheatmap(matrix,main = 'Train data gene expression',
                   fontfamily = "Times",
                   fontsize = 10,
                   annotation_colors = ann.color,
                   annotation_col = a ,
                   annotation_names_col = F,
                   annotation_names_row = F,
                   color = colorRampPalette(c("navy", "white","firebrick1"))(500),
                   cluster_rows = F,
                   cluster_cols = F,
                   show_colnames = F,
                   width = 6,height = 4,
                   file = '05.heatmap.png')

