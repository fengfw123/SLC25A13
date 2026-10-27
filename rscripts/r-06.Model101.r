rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "06_Model101"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)

library(magrittr)
library(tidyr)

# 代码简介：10种具有变量筛选和预后模型构建的算法组合
# 
# 算法包括：Lasso、Ridge、Enet、StepCox、survivalSVM、CoxBoost、SuperPC、plsRcox、RSF、GBM
# 
# 具体内容：
# 在交叉验证的框架下使用一种算法进行变量选择而使用另一种算法构建预后模型，
# 并计算所使用的模型组合(共101种)在外部数据集(或包括训练集)上的一致性指数(C-index)，最终通过热图将模型的评估结果可视化
# 
# 算法输入：
# 训练集表达谱：行为特征(如SYMBOL/ENSEMBL等，但需与测试集有一定交集)，列为样本的表达矩阵，格式见InputData文件中Train_expr.txt
# 训练集生存信息：行为样本，列包括生存状态以及生存时间(代码中需分别指定相对应的变量名)，格式见InputData文件中Train_surv.txt
# 测试集表达谱：行为特征(如SYMBOL/ENSEMBL等，但需与训练集有一定交集)，列为样本的表达矩阵，格式见InputData文件中Test_expr.txt
# 测试集生存信息：行为样本，列包括生存状态以及生存时间，以及一列用于指定队列信息的变量，格式见InputData文件中Test_surv.txt
# 注意：建议对表达谱进行预筛选，保留几十至几千量级的感兴趣特征用于最优模型筛选，否则运行时间将大大延长
# 
# 本代码适用于：
# 存在较多预后意义不明的变量，需同时进行变量筛选，模型筛选和预后模型构建的场景
# 
# 代码输出结果包括：
# 101种算法在训练集上获得的模型(model.rds)，在测试集(或包括训练集)上的评估结果(Cindex_mat)，以及评估结果的热图展示(Cindex.pdf)。
# 
# 参考文献：
# Machine learning-based integration develops an immune-derived lncRNA signature for improving outcomes in colorectal cancer
# DOI: s41467-022-28421-6
# 
# 作者：大鱼海棠
# 工作单位：中国药科大学国家天然药物重点实验室，生物统计与计算药学研究中心
# 目前地址：法国斯特拉斯堡遗传与分子生物研究所（IGBMC），癌症功能基因组研究中心
# 联系邮箱：xlu.cpu@foxmail.com

# 设置工作路径
work.path <- "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/06_Model101/"; setwd(work.path) 

# 设置其他路径
code.path <- file.path(work.path, "Codes") # 存放脚本
raw.path <- file.path(work.path, "rawdata")
data.path <- file.path(work.path, "InputData") # 存在输入数据（需用户修改）
res.path <- file.path(work.path, "Results") # 存放输出结果
fig.path <- file.path(work.path, "Figures") # 存放输出图片

# 如不存在这些路径则创建路径
if (!dir.exists(raw.path)) dir.create(raw.path)
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)


###--- load package ----###
# install.packages("snowfall")
# install.packages("superpc")
# install.packages("survivalsvm")
library(readxl); library(seqinr); library(plyr); library(survival); library(randomForestSRC); library(glmnet);
library(plsRcox); library(superpc); library(gbm); library(mixOmics); library(survcomp); library(CoxBoost);
library(survivalsvm); library(BART); library(snowfall); library(ComplexHeatmap); library(RColorBrewer); library(readr);
library(tidyverse); library(lance); library(tibble); library(openxlsx)


###---- load model file ----####
# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))

# 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiCox表示抽取其他模型所用到的变量并建立多变量cox模型
FinalModel <- c("panML", "multiCox")[2]

## Training Cohort     ---------------------------------------------------------
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，表达谱需有一定变异性，以免建模过程报错)
# load hubgene 
hubgene <- read.csv(file = "../07_Unicox/01.univariate_cox_result_0.05.csv",row.names = 1) %>% rownames()

# train survival data 
# tip: os.time need > zero days 
train.survival <- read.csv(file = "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/SKCM/04.survival.csv")

# train raw data 
train.data <- read.csv(file = "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/SKCM/01.dat.fpkm.csv", header=TRUE, row.names=1, check.names = FALSE)
colnames(train.data) <- gsub('.', '-', colnames(train.data), fixed = T)
train.data <- train.data[, train.survival$sample]

# organize train data 
Train_expr <- train.data[hubgene,] %>% t %>% data.frame()
Train_expr <- scaleData(data = Train_expr, centerFlags = T, scaleFlags = T)

Train_expr <- merge(train.survival, Train_expr, by.x = "sample", by.y = "row.names")
rownames(Train_expr) <- Train_expr$sample
Train_expr <- Train_expr[,-1]

###---- Validation Cohort ----####
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型）
# load validation survival  
survival.va <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/GSE19234/02.survival(GSE19234).csv")

# load validation data 
dat.va <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/GSE19234/01.dat(GSE19234).csv", row.names = 1)
dat.va <- dat.va[, survival.va$sample]

# add cohort 
survival.va$Cohort <- "Testing"

# organize test data 
Test_expr <- dat.va[hubgene,] %>% t %>% data.frame()
Test_expr <- scaleData(data = Test_expr, centerFlags = T, scaleFlags = T) 

Test_expr <- merge(survival.va, Test_expr, by.x = "sample", by.y = "row.names")
rownames(Test_expr) <- Test_expr$sample
Test_expr <- Test_expr[,-1]

# same genes
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))

# 提取相同基因
# identical(colnames(Train_expr),colnames(Test_expr))

# Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
# Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因

# 按队列对数据分别进行标准化（根据情况调整centerFlags和scaleFlags）
## data: 需要表达谱数据（行为样本，列为基因） 
## cohort：样本所属队列，为向量，不输入值时默认全表达矩阵来自同一队列
## centerFlag/scaleFlags：是否将基因均值/标准差标准化为1；
##        默认参数为NULL，表示不进行标准化；
##        为T/F时，表示对所有队列都进行/不进行标准化
##        输入由T/F组成的向量时，按顺序对队列进行处理，向量长度应与队列数一样
##        如centerFlags = c(F, F, F, T, T)，表示对第4、5个队列进行标准化，此时flag顺序应当与队列顺序一致
##        如centerFlags = c("A" = F, "C" = T, "B" = F)，表示对队列C进行标准化，此时不要求flag顺序与data一致
# Train_set <- scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
# # names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) # 注意测试集标准化顺序与此一致
# Test_set <- scaleData(data = Test_expr, centerFlags = T, scaleFlags = T)
# 
# colnames(Train_surv) <- c("OS","OS.time")
# # Test_surv <- Test_surv[,-1]
# colnames(Test_surv) <- c("OS","OS.time","Cohort")
# Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)
# summary(apply(Train_set, 2, var))
# summary(apply(Test_set, 2, var))
# lapply(split(as.data.frame(Test_set), Test_surv$Cohort), function(x) summary(apply(x, 2, var))) # 测试scale结果

###---- Model training and validation ----####
###---- load model list ----####
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有StepCox和RunEnet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)
methods
###---- train model ----#####
# the least cutoff of variable
min.selected.var <- 2 
# define the outcome events that need to be considered, which must appear in the "Train_surv" as well as in the "Test_surv"
timeVar <- "OS.time"; statusVar <- "OS" 

# Pre-training 
Variable <- colnames(Train_expr)[3:length(colnames(Train_expr))]
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))
preTrain.method

# create pre-model
set.seed(seed = 888) # 设置建模种子，使得结果可重复
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # method of machine learning
                                 Train_expr = as.matrix(Train_expr[,3:length(colnames(Train_expr))]), # gene expr of train data 
                                 Train_surv = Train_expr[,1:2], # os and os.time 
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) #  用于训练的生存变量，必须出现在Train_surv中
}
preTrain.var[["simple"]] <- colnames(Train_expr)[3:length(colnames(Train_expr))]

set.seed(seed = 1) # 设置建模种子，使得结果可重复
model <- list() # 初始化模型结果列表
# save.image("preTrain.Rdata")
for (method in methods){ # 循环每一种方法组合
  tryCatch({
    # if error give method debug
    cat(match(method, methods), ":", method, "\n") # 输出当前方法
    method_name <- method # 本轮算法名称
    method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
    
    if (length(method) == 1) method <- c("simple", method)
    
    selected.var = preTrain.var[[method[1]]]
    # 如果筛选出的变量小于阈值，则该算法组合无意义，置空（尤其针对以RSF筛选变量的情况，需在ML脚本中尝试调参）
    if (length(selected.var) <= min.selected.var) {
      model[[method_name]] <- NULL
    } else {
      model[[method_name]] <- RunML(method = method[2], # 用于构建最终模型的机器学习方法
                                    Train_expr = as.matrix(Train_expr[, selected.var]), # 训练集有潜在预测价值的变量
                                    Train_surv = Train_expr[,1:2], # 训练集生存数据
                                    mode = "Model",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                    classVar = classVar)  # 用于训练的生存变量，必须出现在Train_surv中
    }
    
    # 如果最终筛选出的变量小于阈值，则该算法组合也无意义，置空
    if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
      model[[method_name]] <- NULL
    }
  }, error = function(e){
    print('error')
  },finally = {
    print(method)
  })
}

saveRDS(model, file.path(res.path, "model.rds")) # 保存所有模型输出

# # 当要求最终模型为多变量cox时，对模型进行更新
# if (FinalModel == "multiCox"){
#   coxmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量cox模型，从而以cox回归系数和特征表达计算单样本风险得分
#     tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
#                  data = as.data.frame(Train_set[, ExtractVar(fit)]))
#     tmp$subFeature <- ExtractVar(fit) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
#     return(tmp)
#   })
# }
# 
# saveRDS(coxmodel, file.path(res.path, "coxmodel.rds")) # 保存最终以多变量cox拟合所筛选变量的模型

save.image("model.Rdata")

###---- Evaluate the model ----####
# 读取已保存的模型列表（请根据需要调整）
model <- readRDS(file.path(res.path, "model.rds")) # 若希望使用各自模型的线性组合函数计算得分，请运行此行
# model <- readRDS(file.path(res.path, "coxmodel.rds")) # 若希望使用多变量cox模型计算得分，请运行此行

# extract model (variable > cutoff)
methodsValid <- names(model)

# 根据给定表达量计算样本风险评分
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = Test_expr,
                                    type = "lp") # 同原文，使用linear Predictor计算得分
  
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.csv(RS_mat, file.path(res.path, "01.RS_mat.csv"), row.names = T, quote = F) # 输出风险评分文件

RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = Train_expr,
                                    type = "lp") # 同原文，使用linear Predictor计算得分
  
}
RS_train <- as.data.frame(t(do.call(rbind, RS_list)))
write.csv(RS_train, file.path(res.path, "02.RS_train.csv"), row.names = T, quote = F) # 输出风险评分文件


# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]]) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  # 数据框有两列，包含算法以及算法所筛选出的变量
write.csv(fea_df, file.path(res.path, "03.fea_df.csv"), row.names = F, quote = F)

# 对各模型计算C-index
class(Train_expr$OS.time)
class(Test_expr$OS.time)
# Test_surv$OS.time <- as.numeric(Test_surv$OS.time)

Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = as.matrix(Test_expr[,4:length(colnames(Test_expr))]), # 测试集表达谱
                                  Test_surv = Test_expr[,1:3], # 验证集生存数据
                                  Train_expr = as.matrix(Train_expr[,3:length(colnames(Train_expr))]), # 训练集表达谱
                                  Train_surv = Train_expr[,1:2], # 训练集生存数据
                                  cohortVar = "Cohort", # 该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  timeVar = timeVar, #  OS.time
                                  statusVar = statusVar) # OS
}
Cindex_mat <- do.call(rbind, Cindexlist)

write.csv(Cindex_mat, file.path(res.path, "04.cindex_mat.csv"), row.names = T, quote = F)

# Plot --------------------------------------------------------------------
Cindex_mat <- read.csv(file.path(res.path, "04.cindex_mat.csv"), row.names = 1, header = T, check.names = F, stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序

fea_dat <- table(fea_df$algorithm) %>% as.data.frame()
rownames(fea_dat) <- fea_dat$Var1
fea_dat <- fea_dat[names(avg_Cindex), ]
colnames(fea_dat) <- c("Model", "featurenum")
write.csv(fea_dat, file.path(res.path, "05.fea_num.csv"), row.names = T, quote = F)

avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征

# CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置绘图时的队列颜色
CohortCol <- c('seagreen', 'firebrick')
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
cellwidth <- 1; cellheight <- 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "01.heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 7, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())
png(file.path(fig.path, "01.heatmap of cindex.png"), width = cellwidth * ncol(Cindex_mat) + 7, height = cellheight * nrow(Cindex_mat) * 0.45,units = 'in',res = 600)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())


Cindex_mat <- read.table(file.path(res.path, "04.cindex_mat.csv"), row.names = 1, header = T,check.names = F, stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
Cindex_mat <-Cindex_mat [c(1:10),]
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T)
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数

# CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置绘图时的队列颜色
CohortCol <- c('seagreen','firebrick')
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
cellwidth <- 1; cellheight <- 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是s否对行列进行聚类

pdf(file.path(fig.path, "02.heatmap of cindex_top10.pdf"), width = cellwidth * ncol(Cindex_mat) + 7, height = 5)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())
png(file.path(fig.path, "02.heatmap of cindex_top10.png"), width = cellwidth * ncol(Cindex_mat) + 7, height = 5,units = 'in',res = 600)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())