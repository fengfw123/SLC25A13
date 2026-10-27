rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "09_independent"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)

###---- load packages ----####
library(TCGAbiolinks); library(readr); library(lance); library(survival); library(survminer)
library(tidyr); library(rms); library(regplot); library(forestplot); library(cowplot)
###---- load rawdata ----####
risk <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/07_Riskmodel/01.train_riskScore.csv", row.names = 1)
survival <- risk
phenotype <- read_tsv(file = '/data/nas1/yuanshilin/Database/TCGA/TCGA_phenotype/TCGA-SKCM.GDC_phenotype.tsv.gz')
colnames(phenotype)

# train phenotype data 
train.phenotype <- data.frame(id = phenotype$submitter_id.samples,
                              age = phenotype$age_at_initial_pathologic_diagnosis,
                              gender = phenotype$gender.demographic,
                              stage = phenotype$tumor_stage.diagnoses,
                              T.stage = phenotype$pathologic_T,
                              N.stage = phenotype$pathologic_N,
                              M.stage = phenotype$pathologic_M)


# organize data 
train.phenotype <- merge(train.phenotype, survival, by = 'id')
dat <- train.phenotype
rownames(dat) <- dat$id
dat <- dat[,-1]

# age 
table(dat$age, useNA = "ifany")
dat$age <- ifelse(dat$age > 60, ">60", "<=60")

# T stage 
table(dat$T.stage, useNA = "ifany")
dat$T.stage <- gsub('TX',NA,dat$T.stage)
dat$T.stage <- gsub('Tis',NA,dat$T.stage)
dat$T.stage <- gsub("[a-z]", "", dat$T.stage)
# dat$T.stage <- ifelse(dat$T.stage %in% c("T0", "T1"), "T0/1", dat$T.stage)

# N stage 
table(dat$N.stage, useNA = "ifany")
dat$N.stage <- gsub('NX',NA,dat$N.stage,fixed = T)
dat$N.stage <- gsub("[a-z]", "", dat$N.stage)

# M stage 
table(dat$M.stage, useNA = "ifany")
dat$M.stage <- sub("[a-z]",'',dat$M.stage)

# stage 
table(dat$stage, useNA = "ifany")
dat$stage <- gsub("i[^i].*$", "i", dat$stage)
dat$stage <- ifelse(dat$stage == "i" | dat$stage == "not reported", NA, dat$stage)
dat$stage <- ifelse(dat$stage == "stage 0" | dat$stage == "stage i", "stage i/0", dat$stage)

# survival data 
dat$OS.time <- as.numeric(dat$OS.time)
dat$OS <- as.numeric(dat$OS)
dat$group <- factor(dat$risk,levels = c('Low', "High"))

# group 
dat$group <- ifelse(dat$group == "High", "High Risk", "Low Risk")

write.csv(dat, file = '01.phenotype.csv', row.names = F, quote = F)


###---- univariate cox ----####
# factor transformed into character
str(dat)
dat <- lapply(dat, function(x) if(is.character(x)) as.factor(x) else(x))
phenotype.clear <- dat %>% as.data.frame()
str(phenotype.clear)
phenotype.clear$risk <- factor(phenotype.clear$risk,levels = c('Low', "High"))

# univariate cox variates 
coxvariates <- colnames(phenotype.clear)[!colnames(phenotype.clear) %in% c("id", "OS", "OS.time", "riskScore", "group")]

# univariat cox formulas
uncox.formula <- sapply(coxvariates, function(x) {as.formula(paste("Surv(OS.time, OS)~", x))})

# univariate cox model
uncox.model <- lapply(uncox.formula, function(x) {coxph(x, data = phenotype.clear)})


# # # univariate cox ph test
# uncox.zph <- lapply(uncox.model,
#                     function(x) {
#                       y <- cox.zph(x)
#                       #获取p值
#                       p <- y$table[-nrow(y$table),][3]
#                       chisq <- y$table[-nrow(y$table),][1]
#                       df <- y$table[-nrow(y$table),][2]
#                       res <- c(chisq,df,p)
#                       names(res) <- c("chisq","df","p")
#                       return(res)
#                     })
# zph.res <- docall(rbind, uncox.zph) %>% as.data.frame()
# 
# # univariate cox ph plot
# ph.plot <- lapply(uncox.model,
#                   function(model){
#                     zph.test <- cox.zph(model)
#                     zph.plot <- ggcoxzph(zph.test)
#                     zph.plot <- zph.plot[[1]]
#                     return(zph.plot)
#                   })
# 
# combined.plot <- plot_grid(plotlist = ph.plot, ncol = 4)
# ggsave(filename = "02.uncox_ph.pdf", combined.plot, width = 22, height = 10)
# ggsave(filename = "02.uncox_ph.png", combined.plot, width = 22, height = 10, bg = "white", units = "in")

# univariate summary
res.list <- list()
res.ref <- c(1, 1, 1, NA)
for (i in 1:length(uncox.model)) {
  if(length(uncox.model[[i]]$coefficients) > 1){
    tmp.dat <- summary(uncox.model[[i]])
    res.list[[i]] <- rbind(res.ref, cbind(tmp.dat$conf.int[,-2], tmp.dat$coefficients[,5]))
  }else{
    tmp.dat <- summary(uncox.model[[i]])
    res.list[[i]] <- c(tmp.dat$conf.int[-2], tmp.dat$coefficients[5])
  }
}
names(res.list) <- names(uncox.model)
res <- do.call(rbind, res.list) %>% as.data.frame()

# rename
res$Indicators <- c('age', 'gender', 
                   "stage i/0(Reference)", "stage ii","stage iii",
                   'T0(Reference)','T1','T2','T3','T4',
                   # 'T0/1(Reference)','T2','T3','T4',
                   'N0(Reference)','N1','N2','N3',
                   'M1 vs M0', 'risk')
colnames(res) <- c("hr", "low", "up", "pv", "Indicator")

# extract univariate data
res.dat <- data.frame(p.value = res$pv,
                      HR = res$hr,
                      HR.95L = res$low,
                      HR.95H = res$up,
                      row.names = res$Indicator)

write.csv(res.dat, file = "02.univariatecox.csv", quote = F, row.names = T)

###---- forest plot of uncox result ----####
# HR(HR.95L-HR.95H)
hz.dat <- paste(sprintf("%0.2f", res.dat$HR),
                "(",
                sprintf("%0.2f", res.dat$HR.95L),
                "-",
                sprintf("%0.2f", res.dat$HR.95H),
                ")",
                sep = "") 
hz.dat
hz.dat <- ifelse(hz.dat == "1(1-1)", "", hz.dat)

# table text of forest plot 
tabletext <- cbind(c(NA, rownames(res.dat)),
                   c("P value",ifelse(res.dat$p.value < 0.05, "< 0.05", round(res.dat$p.value, 4))),
                   c("Hazard Ratio(95% CI)", hz.dat))
nrow(tabletext) + 1


# extract forest plot 
is.summary <- is.na(res.dat$p.value)

forest.plot <- forestplot(labeltext = tabletext, 
                          graph.pos = 4,  #为Pvalue箱线图所在的位置
                          is.summary = c(TRUE,is.summary),
                          col = fpColors(box="#ed556a", lines="#104680", zero = "gray50"),
                          mean = c(NA,res.dat$HR),
                          lower = c(NA,res.dat$HR.95L), #95%置信区间下限
                          upper = c(NA,res.dat$HR.95H), #95%置信区间上限
                          boxsize = 0.2, lwd.ci = 3,   #箱子大小，线的宽度
                          ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
                          zero = 1,lwd.zero = 0.5,    #zero线宽 基准线的位置
                          colgap = unit(5,"mm"),    #列间隙
                          xticks = c(0,1,2,3,4,18.5), #横坐标刻度
                          clip = c(0, 18.5),
                          lwd.xaxis = 2,            #X轴线宽
                          lineheight = unit(1.2,"cm"), #固定行高
                          graphwidth = unit(.6,"npc"), #图在表中的宽度比例
                          cex = 1.2, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                          hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2),
                                            "18" = gpar(col = "black", lty = 1, lwd = 2)),
                          txt_gp=fpTxtGp(label = gpar(cex = 1),
                                         ticks = gpar(cex = 0.8, fontface = "bold"),
                                         xlab = gpar(cex = 1, fontface = "bold"),
                                         title = gpar(cex = 1.25, fontface = "bold")),
                          xlab = "Hazard Ratio",
                          grid = T,
                          title = "Univariate Cox of independent prognosis") # 垂直于x轴的网格线，对应每个刻度
forest.plot

pdf(file = '03.uncoxforest.pdf', width = 12, height = 8, onefile = F, family = 'Times')
print(forest.plot)
dev.off()

png(file = '03.uncoxforest.png', width = 12, height = 6, units = 'in', res = 600, family = 'Times')
print(forest.plot)
dev.off()


###---- multivariate COX regression ----####
# multivariate variates
mul.variates <- coxvariates[!coxvariates %in% c("gender")]

# multivariate model formula
mul.formula <- as.formula(paste("Surv(OS.time, OS)~", paste0(mul.variates, collapse = "+")))

# multivariate model
mul.model <- coxph(mul.formula, data = phenotype.clear)

# multivariate model ph test 
mulcox.zph <- cox.zph(mul.model)

# multivariate model ph fig
mulph.plot <- ggcoxzph(mulcox.zph)
mulph.plot

# save fig 
combined.plot <- plot_grid(plotlist = mulph.plot, ncol = 3)
ggsave("04.mulcoxph.pdf", combined.plot, width = 16, height = 10)
ggsave("04.mulcoxph.png", combined.plot, width = 16, height = 10, bg = "white", units = "in")

# extract mulcox.zph results
mulcox.table <- mulcox.zph$table[-nrow(mulcox.zph$table),]
mulcox.table

# ph test p > 0.05
res.mul <- summary(mul.model)
res.mul <- cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)


# multivariate result data 
res.mul <- rbind(res.mul[1, ], res.ref, res.mul[c(2:3), ], res.ref, res.mul[4:7, ], res.ref, res.mul[8:12,])
rownames(res.mul)
res.mul$Indicators <- c('age',
                       "stage i/0(Reference)", "stage ii","stage iii",
                       'T0(Reference)','T1','T2','T3','T4',
                       'N0(Reference)','N1','N2','N3',
                       'M1 vs M0', "risk")
colnames(res.mul) = c("hr","low","up","pv","Indicator")


# extract mulresult
mul.res.dat <- data.frame(p.value = res.mul$pv,
                          HR = res.mul$hr,
                          HR.95L = res.mul$low,
                          HR.95H = res.mul$up,
                          row.names = res.mul$Indicator)

write.csv(mul.res.dat, file = "05.mulvariatecox.csv", quote = F, row.names = T)                      

###---- forest plot of mulcox result ----####
# HR(HR.95L-HR.95H) 
mulhz.dat <- paste(sprintf("%0.2f",mul.res.dat$HR),
                   "(",
                   sprintf("%0.2f",mul.res.dat$HR.95L),
                   "-",
                   sprintf("%0.2f",mul.res.dat$HR.95H),
                   ")",
                   sep = "")
                
mulhz.dat
mulhz.dat <- ifelse(mulhz.dat == "1(1-1)", "", mulhz.dat)

# table text of forest plot 
mul.tabletext <- cbind(c(NA, rownames(mul.res.dat)),
                       c("P value", ifelse(mul.res.dat$p.value < 0.05, "< 0.05", round(mul.res.dat$p.value, 4))),
                       c("Hazard Ratio(95% CI)", mulhz.dat))
nrow(mul.tabletext) + 1

# extract forest plot 
is.summary <- is.na(mul.res.dat$p.value)

mul.forest.plot <- forestplot(labeltext = mul.tabletext,
                              graph.pos = 4,  #为Pvalue箱线图所在的位置
                              is.summary = c(TRUE,is.summary),
                              col = fpColors(box="#ed556a", lines="#104680", zero = "gray50"),
                              mean =c(NA, mul.res.dat$HR),
                              lower = c(NA, mul.res.dat$HR.95L), #95%置信区间下限
                              upper = c(NA, mul.res.dat$HR.95H), #95%置信区间上限
                              boxsize = 0.2, lwd.ci = 3,   #箱子大小，线的宽度
                              ci.vertices.height = 0.08, ci.vertices=TRUE, #置信区间用线宽、高、型
                              zero = 1,lwd.zero = 0.5,    #zero线宽 基准线的位置
                              colgap = unit(5,"mm"),    #列间隙
                              xticks = c(0,1,2,3,4,27), #横坐标刻度
                              clip = c(0, 18.5),
                              lwd.xaxis = 2,            #X轴线宽
                              lineheight = unit(1.2,"cm"), #固定行高
                              graphwidth = unit(.6,"npc"), #图在表中的宽度比例
                              cex = 1.2, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                              hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2),
                                                "17" = gpar(col = "black", lty = 1, lwd = 2)),
                              txt_gp = fpTxtGp(label = gpar(cex = 1),
                                               ticks = gpar(cex = 0.8, fontface = "bold"),
                                               xlab = gpar(cex = 1, fontface = "bold"),
                                               title = gpar(cex = 1.25, fontface = "bold")),
                              xlab = "Hazard Ratio",
                              grid = T,
                              title = "multvariate Cox of independent prognosis") # 垂直于x轴的网格线，对应每个刻度
                          
mul.forest.plot

pdf(file = '06.mulcoxforest.pdf', width = 12, height = 8, onefile = F, family = 'Times')
print(mul.forest.plot)
dev.off()

png(file = '06.mulcoxforest.png', width = 12, height = 6, units = 'in', res = 600, family = 'Times')
print(mul.forest.plot)
dev.off()


# ###---- nomogram ----####
# no.formula <- mul.formula
# 
# # create nomogram
# ddist <- datadist(phenotype.clear)
# options(datadist = 'ddist')
# res.no <- psm(no.formula, data = phenotype.clear, dist = 'lognormal')
# 
# # survival function 
# surv <- Survival(res.no)
# 
# #function(x) surv(365, x) # 1年事件发生概率
# #function(x) surv(730, x) # 2年事件发生概率
# function(x) surv(1095, x) # 3年事件发生概率
# function(x) surv(1825, x) # 5年事件发生概率
# function(x) surv(2555, x) # 7年事件发生概率
# nom.cox <- nomogram(res.no,
#                     fun = list(function(x) surv(1095, x),
#                                function(x) surv(1825, x),
#                                function(x) surv(2555, x)),
#                     funlabel=c("3-year Survival Probability", "5-year Survival Probability", "7-year Survival Probability"),
#                     maxscale = 10,
#                     fun.at = c(0.01, seq(0.1, 0.9, by = 0.2), 0.95, 0.99),
#                     lp = F)
# nocox <- cph(no.formula, data = phenotype.clear, x = TRUE, y = TRUE, surv = TRUE, na.action = na.delete)
# 
# 
# regplot(nocox, plots = c("violin", "bars"),
#              observation = T, title = "Survival Nomogram",
#              failtime = c(365*3, 365*5,365*7), prfail = T,  points = TRUE)
# save.p <- recordPlot()
# 
# pdf("06.nomogram.pdf", width = 14, height = 10, family = 'Times')
# replayPlot(save.p)
# dev.off()
# 
# png("06.nomogram.png", width = 14, height = 10, family = 'Times', units = 'in', res = 600)
# replayPlot(save.p)
# dev.off()









