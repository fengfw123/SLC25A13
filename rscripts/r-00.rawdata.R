rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "00_rawdata"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)


pacman::p_load(GEOquery, Biobase, readr, readxl, tidyverse, dplyr)
library(AnnoProbe) #idmap
library(TCGAbiolinks)
################################# 数据集准备 ###################################
### 训练集 TCGA-SKCM 预后
### 训练集 GSE15605  候选基因（差异）
### 验证集 GSE65904  预后验证 GSE19234
### 单细胞 GSE222446 单细胞
### 线粒体 MDGs      方向1
### BRAF   BRAFGs    方向2

# 01. TCGA-SKCM SKCM  ----------------------------------------------------------
gse <- "SKCM"
setwd("/data/nas1/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/")
if(!dir.exists(gse)){
  dir.create(gse)
}
setwd(gse)

## 读取从xena下载的数据
tcga.expr <- read_tsv(file = '/data/nas1/yuanshilin/Database/TCGA/TCGA.matrix/TCGA-SKCM.htseq_counts.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr <- 2^tcga.expr-1

## 对数据进行id转化
genecode <- read.table(file = '/data/nas1/yuanshilin/Database/TCGA/Gencode/gencode.v22.annotation.gene.probeMap')
probe2symbol <- genecode[,(1:2)]
colnames(probe2symbol) <- c('ID','symbol')
probe2symbol <- probe2symbol[-1,]

dat.tcga <- tcga.expr
dat.tcga$ID <- rownames(dat.tcga) 
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.tcga)
### 58387 472
## 筛选癌症组织，去掉癌旁组织。01-09为肿瘤，10-19为正常对照
mete <-data.frame(colnames(dat.tcga))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}

names(mete)=c("id","group")
table(mete$group)
### Control 1 Trum 471

mete$group <- as.factor(mete$group)
mete <- subset(mete,mete$group=="T")

## 保留有生存数据的样例
exp_tumor<-dat.tcga[,which(colnames(dat.tcga)%in%mete$id)]
survival <- read.delim2('/data/nas1/yuanshilin/Database/TCGA/TCGA_survival/TCGA-SKCM.survival.tsv')
survival <- survival[,-3]
exp_tumor <- exp_tumor[,colnames(exp_tumor)%in%survival$sample]
survival <- survival[survival$sample %in% colnames(exp_tumor),]
dim(exp_tumor)
### 58387   457个

### 预后保留编码蛋白得基因
pcg <- read.delim2('/data/nas1/yuanshilin/Database/TCGA/PCG.xls(v22)')
dat.final <- exp_tumor[pcg$gene_name,]
dat.final <- na.omit(dat.final)
dim(dat.final)
### 19814   457个

##fpkm
expr_fpkm <- read_tsv(file = '/data/nas1/yuanshilin/Database/TCGA/TCGA.matrix/TCGA-SKCM.htseq_fpkm.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
## xena下载的数据经过了log2+1转化，需要将其还原
expr_fpkm<-2^expr_fpkm-1
## 对数据进行id转化
dat_fpkm <- expr_fpkm
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID <- as.character(dat_fpkm$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)
dat_fpkm <- dat_fpkm %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])  ## 把第一列变成行名并删除
colnames(dat_fpkm)[1:5]
# colnames(dat_fpkm) <- gsub('.', '-', colnames(dat_fpkm), fixed = T)
dim(dat_fpkm)


dat_fpkm <- dat_fpkm[,colnames(dat.final)]
dat_fpkm <- dat_fpkm[rownames(dat.final),]
dat_fpkm <- na.omit(dat_fpkm)
dim(dat_fpkm)
### 457


## 输出fpkm文件
write.csv(dat_fpkm, file = '01.dat.fpkm.csv', row.names = T, quote = F)
## 输出counts文件
dat.final <- dat.final[rownames(dat_fpkm),]
dim(dat.final)

write.csv(dat.final, file = '02.dat.tcga.csv', row.names = T, quote = F)


## fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
dim(dat_tpm)
dat_tpm <- as.data.frame(dat_tpm)
dat_tpm <- log2(dat_tpm + 1)
write.csv(dat_tpm, file = '03.dat.tpm(log).csv', row.names = T, quote = F)

write.csv(survival, file = "04.survival.csv", row.names = F)

#临床特征  ---------------------------------------------------------------------
phe <- read.delim2("/data/nas1/yuanshilin/Database/TCGA/TCGA_phenotype/TCGA-SKCM.GDC_phenotype.tsv.gz", sep = "\t")
cli_dat <- data.frame(sample = phe$submitter_id.samples,
                      gender = phe$gender.demographic,
                      age = phe$age_at_initial_pathologic_diagnosis,
                      stage = phe$tumor_stage.diagnoses,
                      T_stage = phe$pathologic_T,
                      N_stage = phe1$pathologic_N,
                      M_stage = phe1$pathologic_M
)
cli_dat[cli_dat == ""] <- NA
cli_dat <- na.omit(cli_dat)
#cli_dat$age <- as.numeric(cli_dat$age) / 365
write.csv(cli_dat, file = "05.phenotype.csv", quote = F, row.names = F)


# 02. GSE15605 -----------------------------------------------------------------
rm(list = ls())
gse <- "GSE15605"
setwd("/data/nas1/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/")
if(!dir.exists(gse)){
  dir.create(gse)
}
setwd(gse)

### 获取数据压缩包
gset <- getGEO(gse,
               destdir = '.',
               GSEMatrix = T,
               getGPL = F)  #为了下载速度不下载平台文件 

### 提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a <- gset[[1]]

### 提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
expr <- as.data.frame(exprs(gset[[1]])) 

### 判断下是否需要log
exprSet <- expr
ex <- exprSet #把需要判断的矩阵命名为ex
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

gset$GSE15605_series_matrix.txt.gz@annotation
### GPL570

### 提取芯片平台信息
gpl <- getGEO("GPL570",destdir = '/data/nas1/yuanshilin/Pipline/GPL/')
gpl <- Table(gpl)
colnames(gpl)

### 提取基因Symbol和ID
probe2symobl <- gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
### 根据数据调整
probe2symobl <- probe2symobl[probe2symobl$symbol!='',]
# probe2symobl <- subset(probe2symobl, probe2symobl$symbol != "" & probe2symobl$symbol != "---")

### 与表达矩阵匹配
dat <- expr
dat$ID <- rownames(dat)
dat$ID <- as.character(dat$ID)
probe2symobl$ID <- as.character(probe2symobl$ID)

### 整理数据
dat <- dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

### 提取样本信息
pd <- pData(a)
colnames(pd)
table(pd$`diagnosis:ch1`)

### 提取分组信息矩阵
group <- data.frame(sample = pd$geo_accession, group = pd$source_name_ch1,
                    subgroup = pd$`diagnosis:ch1`)
table(group$subgroup)
group <- group[group$subgroup %in% c("NORMAL SKIN", "SSM"), ] %>% select(sample, group)
table(group$group)
# Frozen half of the melanoma punch biopsy
# 28                

### 选择“疾病样本”和“正常样本”
group$group <- ifelse(group$group=='Frozen half of the melanoma punch biopsy','SKCM','Control')
table(group$group)
# Control    SKCM 
# 16      28

### 排序,将Control or Normal放在上边 
group <- group[order(group$group, decreasing = F),]

### 筛选所需样本的表达矩阵
dat <- dat[,group$sample]
dat$gene <- rownames(dat)
dat2 <- dat

### 去除空格和重复
dat2$gene <- gsub(' ','',dat2$gene)
dat2 <- dat2[!duplicated(dat2$gene),]

### 更改行名，去除ID列
rownames(dat2) <- dat2$gene
dat2 <- dat2[,which(colnames(dat2) != "gene")]

write.csv(dat2, file = paste0('01.dat(', gse ,').csv'), quote = F, row.names = T)
write.csv(group, file = paste0('02.group(', gse ,').csv'), quote = F, row.names = F)



# 02-2.GSE19234 ------------------------------------------------------------
rm(list = ls())
gse <- "GSE19234"
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/")
if(!dir.exists(gse)){
  dir.create(gse)
}
setwd(gse)
gset <- getGEO("GSE19234",
               destdir = '.',
               GSEMatrix = T,
               getGPL = F)  #为了下载速度不下载平台文件 
# class(gset)

#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a <- gset[[1]] 
# class(a)

#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat <- as.data.frame(exprs(gset[[1]])) 

gset$GSE19234_series_matrix.txt.gz@annotation
### GPL570
gpl <- getGEO("GPL570",destdir = '/data/nas1/yuanshilin/project/01.CSZK20208-1/00_rawdata/')
gpl <- Table(gpl)    

colnames(gpl)
probe2symobl <- gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='') %>% 
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl <- probe2symobl[probe2symobl$symbol!='',]

dat$ID <- rownames(dat)
dat$ID <- as.character(dat$ID)
probe2symobl$ID <- as.character(probe2symobl$ID)
dat <- dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

#判断且进行log转化
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex + 1)
print("log2 transform finished")}else{print("log2 transform not needed")}

dat <- ex 
#df <- log2(df+1) 
# write.csv(dat,file = paste0('01.dat(GSE19234.csv'),quote = F)
pd <- pData(a) #查看a的临床信息，为后面选择用于分组的变量做准备
colnames(pd)

#survive<- pd[,grepl("days since initial diagnosis:ch1|staus dead or alive:ch1", colnames(pd))]  #提取指定的列名
#setwd("../10_External_va")
#write.table(survive,file = 'survive.csv',sep = ',',row.names = T,quote = F)
survival.va <- data.frame(sample=pd$geo_accession,OS=pd$`staus dead or alive:ch1`,OS.time=pd$`days since initial diagnosis:ch1`)
write.csv(dat, file = paste0('01.dat(', gse ,').csv'), quote = F, row.names = T)
write.csv(survival.va, file = paste0('02.survival(', gse ,').csv'), quote = F, row.names = F)


# 03. GSE65904  ----------------------------------------------------------------
rm(list = ls())
gse <- "GSE65904"
setwd("/data/nas1/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/")
if(!dir.exists(gse)){
  dir.create(gse)
}
setwd(gse)


### 获取数据压缩包
gset <- getGEO(gse,
               destdir = '.',
               GSEMatrix = T,
               getGPL = F)  #为了下载速度不下载平台文件 

### 提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a <- gset[[1]]

### 提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
expr <- as.data.frame(exprs(gset[[1]])) 

### 判断下是否需要log
exprSet <- expr
ex <- exprSet #把需要判断的矩阵命名为ex
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")
expr <- exprSet}else{print("log2 transform not needed")}

gset$GSE65904_series_matrix.txt.gz@annotation
### GPL10558

### 提取芯片平台信息
gpl <- getGEO(gset$GSE65904_series_matrix.txt.gz@annotation,
              destdir = '/data/nas1/yuanshilin/Pipline/GPL/')
# gpl <- Table(gpl)
# colnames(gpl)
gpl <- idmap(gpl = 'GPL10558')
### 提取基因Symbol和ID
probe2symobl <- gpl
# dplyr::select('ID','Symbol')%>%
# filter('Symbol'!='')%>%
# separate('Symbol',c('symbol','drop'),sep = '///')%>%
# dplyr::select(-drop)
### 根据数据调整
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
colnames(probe2symobl)<-c('ID','symbol')

# probe2symobl <- probe2symobl[probe2symobl$symbol!='',]
# probe2symobl <- subset(probe2symobl, probe2symobl$symbol != "" & probe2symobl$symbol != "---")

### 与表达矩阵匹配
dat <- expr
dat$ID <- rownames(dat)
dat$ID <- as.character(dat$ID)
probe2symobl$ID <- as.character(probe2symobl$ID)

### 整理数据
dat <- dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

### 提取样本信息
pd <- pData(a)
colnames(pd)
table(pd$`tumor stage:ch1`)
### 提取生存信息
survival <- data.frame(sample = pd$geo_accession,
                       OS = pd$`disease specific survival (1=death, 0=alive):ch1`,
                       OS.time = pd$`disease specific survival in days:ch1`)

### 提取SKCM样本
group <- data.frame(sample = pd$geo_accession,
                    group = pd$`tissue:ch1`)

# group$group <- substr(group$group, 1, 17)
table(group$group)
# Cutaneous
#  22
### 选择皮肤样本
group <- group[group$group == "Cutaneous",]

### 筛选有生存数据的表达矩阵
survival <- survival[survival$sample %in% group$sample,]
survival <- survival[!survival$OS == "NA",]
dat <- dat[,survival$sample]
dat$gene <- rownames(dat)
dat2 <- dat

### 去除空格和重复
dat2$gene <- gsub(' ','',dat2$gene)
dat2 <- dat2[!duplicated(dat2$gene),]

### 更改行名，去除ID列
rownames(dat2) <- dat2$gene
dat2 <- dat2[,which(colnames(dat2) != "gene")]
# dat2 <- log2(dat2+6)

write.csv(dat2, file = paste0('01.dat(', gse ,').csv'), quote = F, row.names = T)
write.csv(survival, file = paste0('02.survival(', gse ,').csv'), quote = F, row.names = F)


# 04. MDGs  --------------------------------------------------------------------
rm(list = ls())
gse <- "MDGs"
setwd("/data/nas1/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/")
if(!dir.exists(gse)){
  dir.create(gse)
}
setwd(gse)
MDGs <- data.frame(symbol = c("ARMC10", "DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2", "MTFP1", "MTFR1",
                              "MTFR2", "MUL1", "OMA1", "RAB24", "SLC25A46", "SPIRE1", "STX17", "ARL2", 
                              "MFN1", "MFN2", "MIGA1", "MIGA2", "MTCH2", "OPA1", "PLD6"))
write.csv(MDGs, file = "01.Mitochondrial.csv", quote = F, row.names = F)


# 05. BRAFGs  ------------------------------------------------------------------
rm(list = ls())
gse <- "BRAFs"
setwd("/data/nas1/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/")
if(!dir.exists(gse)){
  dir.create(gse)
}
setwd(gse)

BRAFs <- data.frame(symbol = c("SNX10","S100B","CSPG4","AIM2","FDCSP","KCNJ10",
                               "COL22A1","ITIH6","IL27RA","SERPINE2","OLIG1",
                               "CD200","METTL7B","TIMP1"))

write.csv(BRAFs, file = "01.BRAFs.csv", quote = F, row.names = F)



