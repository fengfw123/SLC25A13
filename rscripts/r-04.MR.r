rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "04_MR"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)

# fc = 2
# mr分析-------------
if (TRUE){
  if (! dir.exists("./00.steiger_dat.all")){
    dir.create("./00.steiger_dat.all")
  }
  
  if (! dir.exists("./00.result.all")){
    dir.create("./00.result.all")
  }
  
  if (! dir.exists("./00.R&F.all")){
    dir.create("./00.R&F.all")
  }
  
  if (! dir.exists("./00.heterogeneity.all")){
    dir.create("./00.heterogeneity.all")
  }
  
  if (! dir.exists("./00.pleiotropy.all")){
    dir.create("./00.pleiotropy.all")
  }
  
  # if (! dir.exists("./00.presso.all")){
  #   dir.create("./00.presso.all")
  # }
}





library(clusterProfiler)
Candidated_genes<-read.csv('../02_Hubgene/02.hubgenes.csv', row.names = 1)

id <- bitr(Candidated_genes$symbol,fromType = "SYMBOL",
         toType = c("ENSEMBL" ),
         OrgDb = 'org.Hs.eg.db')
id$ENSEMBL <- paste('eqtl','a',id$ENSEMBL,sep='-')
ids <- id$ENSEMBL 

#BiocManager::install("VariantAnnotation")
#BiocManager::install("gwasglue")
#devtools::install_github("huangyebao/ClumpMR")
#devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
library(VariantAnnotation)
library(gwasglue)
library(data.table)
library(tidyr )
library(dplyr)
library(ieugwasr)
library(TwoSampleMR)
library(ClumpMR)
library(ggplot2)
###################### outcome vcf
# # finn文件  从芬兰数据库搜索，https://storage.googleapis.com/finngen-public-data-r9/summary_stats/ + 编号
# outcome_file <- "finngen_R9_C3_MELANOMA_SKIN_EXALLC.gz"
# outcome<-data.table::fread(outcome_file,sep = "\t",data.table=F)
# colnames(outcome) <- c("chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls")
# outcome$pval = 10^(-outcome$pval)
# outcome_all <- format_data(
#   dat=outcome,
#   type = "outcome",
#   snps = NULL,
#   header = TRUE,
#   snp_col = "rsids",
#   beta_col = "beta",
#   se_col = "sebeta",
#   eaf_col = "af_alt",
#   effect_allele_col = "alt",
#   other_allele_col = "ref",
#   pval_col = "pval",
#   # ncase_col = "n_cases",
#   # ncontrol_col = "n_controls",
#   chr_col = "chrom",
#   pos_col = "pos"
# )
# outcome$ncase.exposure = 98
# outcome$ncontrol.exposure = 174006
# save(outcome, file='outcome.RData')


# #vcf文件  从opengwas数据库下载
# vcf2 <- VariantAnnotation::readVcf('ieu-a-32.vcf.gz', "hg19")
# outcome <- gwasvcf_to_TwoSampleMR(vcf2)
# outcome$ncase.exposure = 6968
# outcome$ncontrol.exposure = 20464
# save(outcome, file='outcome.RData')
load('outcome.RData')

# tsv文件 从gwas catelog下载
# outcome<-read.table('./ebi-a-GCST003484.vcf.gz',header=F)     ##ebi-a-GCST003484
# outcome<-VariantAnnotation::readVcf('ebi-a-GCST003484.vcf.gz', "hg19")
# outcome <- gwasvcf_to_TwoSampleMR(outcome)
# save(outcome, file='outcome.RData')
# load('/data/nas1/liky/project/38_BJTC-574-2/07_MR_kb10/outcome.RData')

#outcome <- gwasvcf_to_TwoSampleMR(outcome, type='outcome')
list_gene <- read.table("/data/nas2/database/MR/eQTL_vcf/eQTL_list")
ids <- id$ENSEMBL[id$ENSEMBL %in% list_gene$V1]
#save.image(file="MRdat.RData")
# load('MRdat.RData')
# load('/data/nas1/liky/project/38_BJTC-574-2/07_MR_kb10/outcome.RData')
#results_all<-data.frame()


# i <- 60
for (i in 1:length(ids))  {
  #for (i in c(2,10,13))  {
  tryCatch({
    #{i <- 21
    setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
    dir.create(paste(ids[i],sep='_'))
    setwd(paste(ids[i],sep='_'))
    
    
    ##暴露因素
    vcf1 <- VariantAnnotation::readVcf(paste('/data/nas2/database/MR/eQTL_vcf/',ids[i],'.vcf.gz',sep=''),"hg19")
    # 将vcf转化为TwoSampleMR的格式
    exp_data <- gwasvcf_to_TwoSampleMR(vcf1)
    
    temp_dat <- exp_data[exp_data$pval.exposure < 5e-07,] #暴露因素和SNP强相关
    # 补充ld_clump需要的三列
    temp_dat$id <- temp_dat$id.exposure # id列用于区分是单个或多个来源的gwas数据，clump按照单个gwas数据进行
    temp_dat$rsid <- temp_dat$SNP
    temp_dat$pval <- temp_dat$pval.exposure
    # 执行ld_clump
    # bfile指定参考文件的路径
    exp_dat <- ld_clump(temp_dat,
                        plink_bin = get_plink_exe(),
                        bfile = '/data/nas2/database/MR/Coloc_bfile/EUR' ,
                        clump_kb = 1000, clump_r2 = 0.001)
    exp_dat$id<-ids[i]
    exp_dat$exposure = paste(ids[i],' || id:finngen_R9_C3_MELANOMA_SKIN_EXALLC',sep='')
    #exp_dat <- exp_dat[exp_dat$mr_keep.exposure=='TRUE',]
    write.csv(exp_dat,'01.exposure_dat_input.csv',quote=F,row.names=F)
    # ## 解决变量 table文件
    # outcome_dat <- format_data(
    #   dat=outcome,
    #   type = "outcome",
    #   snps = exp_dat$rsid,
    #   header = TRUE,
    #   snp_col = "variant_id",
    #   beta_col = "beta",
    #   se_col = "standard_error",
    #   effect_allele_col = "effect_allele",
    #   other_allele_col = "other_allele",
    #   pval_col = "p_value",
    #   ncase_col = "n_cases",
    #   ncontrol_col = "n_controls",
    #   chr_col = "chromosome",
    #   pos_col = "base_pair_location"
    # )
    #VCF文件
    outcome_dat <- format_data(
      dat = outcome,
      type = "outcome",
      snps = exp_dat$SNP,
      header = TRUE,
      snp_col = "rsids",
      beta_col = "beta",
      se_col = "sebeta",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      pval_col = "pval",
      ncase_col = "ncase.exposure",
      ncontrol_col = "ncontrol.exposure",
      chr_col = "chrom",
      pos_col = "pos"
    )
    ## ## 解决变量 vcf文件  
    #outcome_dat<- subset(outcome, outcome$SNP%in%exp_dat$SNP)
    
    outcome_dat$id.outcome <- "SKCM"
    outcome_dat$outcome <- "SKCM || finngen_R9_C3_MELANOMA_SKIN_EXALLC"
    write.csv(outcome_dat,'02.outcome_dat_input.csv',quote=F,row.names=F)
    dat <- harmonise_data(exposure_dat =exp_dat,outcome_dat = outcome_dat) #表格里有回文SNP，没有删除，后续分析会根据MR_keep为T自动筛选
    dat <- dat[dat$mr_keep=='TRUE',]
    if(nrow(dat)<3){next}
    write.csv(dat,file = '03.harmonized_dat_input.csv',row.names=F)
    
    
    #steiger_sl<-steiger_filtering(dat) #判断单一SNP的方向性是否正确，SNP与暴露强相关，与结局无关
    steiger_sl <- directionality_test(dat)
    write.csv(steiger_sl,file = '04.steiger_dat.csv',row.names=F)
    write.csv(steiger_sl,file = paste0("../00.steiger_dat.all/",ids[i],'_',i,'_steiger_dat.csv'),row.names=F)
    
    result <- generate_odds_ratios(mr_res = mr(dat))  
    #result <- generate_odds_ratios(mr_res = mr(dat, method_list = c("mr_ivw_mre")))   
    write.csv(result,file = '05.result_method.csv',row.names=F)
    write.csv(result,file = paste0("../00.result.all/",ids[i],'_',i,'_result_method.csv'),row.names=F)
    
    
    clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
    clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
    #clumping <- clumping[clumping$F>10,] #F值 弱工具变量 要大于10
    write.csv(clumping,file = '06.R&Fstat.csv',row.names = F,quote = F)
    write.csv(clumping,file = paste0("../00.R&F.all/",ids[i],'_',i,'_clumping.csv'),row.names=F)
    
    
    library(patchwork)
    library(ggplot2)
    p1 <- mr_scatter_plot(dat,mr_results = mr(dat,method_list = c('mr_ivw')))[[1]] +
      scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.position = 'top',
            legend.direction = 'vertical',
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('01.scatter.pdf',w=12,h=12,family='Times')
    print(p1)
    dev.off()
    
    png('01.scatter.png',w=12,h=12,units='in',res=600,family='Times')
    print(p1)
    dev.off()
    
    #森林图
    res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
    p2<-mr_forest_plot(res_single)[[1]] +
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            # legend.title = element_blank(),
            # legend.position = c(0.1,0.9),
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.title = element_blank(),
            legend.position = 'none',
            legend.direction = 'vertical',
            plot.title=element_text(size=8),
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('02.forest.pdf',w=18,h=20,family='Times')
    print(p2)
    dev.off()
    
    png('02.forest.png',w=18,h=20,units='in',res=600,family='Times')
    print(p2)
    dev.off()
    
    p3<-mr_funnel_plot(res_single)[[1]] +
      scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.position = 'top',
            legend.direction = 'vertical',
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('03.funnel.pdf',w=10,h=11,family='Times')
    print(p3)
    dev.off()
    
    png('03.funnel.png',w=10,h=11,units='in',res=600,family='Times')
    print(p3)
    dev.off()
    
    #异质性检验
    het <- mr_heterogeneity(dat)
    het
    write.csv(het,file = '10.heterogeneity.csv',row.names=F)
    write.csv(het,file = paste0("../00.heterogeneity.all/",ids[i],'_',i,'_heterogeneity.csv'),row.names=F)
    
    
    ##多效性检测
    ple <- mr_pleiotropy_test(dat)
    #<0.05就别做了
    ple
    write.csv(ple,file = '11.pleiotropy.csv',row.names=F)
    write.csv(ple,file = paste0("../00.pleiotropy.all/",ids[i],'_',i,'_pleiotropy.csv'),row.names=F)
    
    
    #逐个剔除检验
    single <- mr_leaveoneout(dat)
    p4<-mr_leaveoneout_plot(single)[[1]] +
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.title = element_blank(),
            legend.position = 'none',
            legend.direction = 'vertical',
            plot.title=element_text(size=8),
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('04.leaveoneout.pdf',w=18,h=18,family='Times')
    print(p4)
    dev.off()
    
    png('04.leaveoneout.png',w=18,h=18,units='in',res=600,family='Times')
    print(p4)
    dev.off()
    
    # library(MRPRESSO) #要满足p>0.05
    # mr_presso<-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
    # presso <- data.frame(id=paste(ids[i],' || id:ieu-a-32',sep=''), Pvalue=mr_presso$`MR-PRESSO results`$`Global Test`$Pvalue)
    # write.csv(presso,file = '13.mr_presso.csv',row.names=F)
    # write.csv(presso,file = paste0("../00.presso.all/",ids[i],'_',i,'_presso.csv'),row.names=F)
    #}
    
    
  }, error = function(e){
    print('error')
  },finally = {
    print(ids[i])
    print(i)
  })
}

#初步整理表格-----------
rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "04_MR"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)

library(clusterProfiler)
Candidated_genes <- read.csv('../02_Hubgene/02.hubgenes.csv',row.names = 1)

id <- bitr(Candidated_genes$symbol,fromType = "SYMBOL",
         toType = c("ENSEMBL"),
         OrgDb = 'org.Hs.eg.db')
id$ENSEMBL <- paste('eqtl','a',id$ENSEMBL,sep='-')
ids <- id$ENSEMBL

path <- getwd()
herb <- list.files(path = "./00.result.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.result.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
results_all <- Reduce(rbind,tem)
results_all$id <- results_all$exposure
results_all <- results_all[,c('id','method','nsnp','pval','or','or_lci95','or_uci95')]
#write.csv(results_all, '00.results_all.csv')
# results_all <- results_all[results_all$method=='Inverse variance weighted' | 
#                              results_all$method=='Inverse variance weighted (multiplicative random effects)',]

setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
path <- getwd()
herb <- list.files(path = "./00.R&F.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.R&F.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
RF_all <- Reduce(rbind,tem)
RF_all$id <- RF_all$exposure
RF_all <- RF_all[,c('id','SNP','mr_keep','F')]
#write.csv(heterogeneity_all, '00.heterogeneity_all.csv')


setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
path <- getwd()
herb <- list.files(path = "./00.heterogeneity.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.heterogeneity.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
heterogeneity_all <- Reduce(rbind,tem)
heterogeneity_all <- heterogeneity_all[heterogeneity_all$method=='Inverse variance weighted',]
heterogeneity_all$id <- heterogeneity_all$exposure
heterogeneity_all <- heterogeneity_all[,c('id','Q_pval')]
colnames(heterogeneity_all) <- c('id','heterogeneity_pval')
#write.csv(heterogeneity_all, '00.heterogeneity_all.csv')

setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
path <- getwd()
herb <- list.files(path = "./00.pleiotropy.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.pleiotropy.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
pleiotropy.all <- Reduce(rbind,tem)
pleiotropy.all$id <- pleiotropy.all$exposure
pleiotropy.all <- pleiotropy.all[,c('id','pval')]
colnames(pleiotropy.all) <- c('id','pleiotropy_pval')

setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
path <- getwd()
herb <- list.files(path = "./00.steiger_dat.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.steiger_dat.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
steiger_dat.all <- Reduce(rbind,tem)
steiger_dat.all$id <- steiger_dat.all$exposure
#steiger_dat.all <- steiger_dat.all[,c('id','SNP','steiger_dir','steiger_pval')]
steiger_dat.all <- steiger_dat.all[,c('id','correct_causal_direction','steiger_pval')]

# setwd("/data/nas1/liky/project/38_BJTC-574-2/09_MR_ieu32_kb10")
# path <- getwd()
# herb <- list.files(path = "./00.presso.all/")
# tem <- lapply(herb,function(x){
#   tw <- read.csv(paste0("./00.presso.all/",x))
#   #tw$herb = gsub('.csv','',x)
#   tw
# })
# presso.all <- Reduce(rbind,tem)
#presso.all$id <- presso.all$exposure
#presso.all <- presso.all[,c('id','SNP','steiger_dir','steiger_pval')]

result.write1 <- merge(results_all, heterogeneity_all, by='id')
result.write1 <- merge(result.write1, pleiotropy.all, by='id')
result.write1 <- merge(result.write1, steiger_dat.all, by='id')
result.write1 <- merge(result.write1, RF_all, by='id')


result.write <- result.write1
result.write$ENSEMBL<-substr(result.write$id,1,22)
result.write<-merge(result.write,id,by='ENSEMBL')
#result.write <- merge(result.write,presso.all,by='id',all.x=TRUE)
write.csv(result.write,'01.results_all.csv')

result.write2<-result.write[result.write$method=='Inverse variance weighted',]
result.write2<-result.write2[result.write2$pval<0.05,]
write.csv(result.write2, "01.results_0.05.csv")
length(result.write2$SYMBOL %>% unique())
# 46
result.write <- read.csv('01.results_all.csv',row.names = 1)
length(result.write$SYMBOL %>% unique())
# 108



#异质性检验基因--------------
rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "04_MR"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)
getwd()


# cox <- read.csv('sig_diff_write_25_GSE16879.csv')
# cox <- cox[,c(1,2)]

mr <- read.csv('01.results_all.csv')
table(mr$mr_keep)
mr[order(mr$F),]$F[1]

a <- mr[,c(3,4,6,7,10,17)] %>% na.omit() %>% unique()
a <- a[a$heterogeneity_pval<0.05,]
a <- a[a$method == "Inverse variance weighted", ]
a$id <- gsub(" \\|\\| id:finngen_R9_C3_MELANOMA_SKIN_EXALLC", "", a$id)


herb <- list.files(path = "./00.result.all/") %>% as.data.frame()
herb <- separate(herb, ".", c("id", "i", "drop2"), sep = "_")[,1:2]
a <- merge(a, herb, by = "id")
a <- a$i %>% as.numeric()
  
  
#异质性检验基因mr_ivw_mre重新跑--------------
# rm(list = ls())
setwd('/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/')
if (! dir.exists("./00.mr_ivw_mre")){
  dir.create("./00.mr_ivw_mre")
}
setwd("./00.mr_ivw_mre")

if (TRUE){
  
  if (! dir.exists("./00.result.all")){
    dir.create("./00.result.all")
  }
  
  if (! dir.exists("./00.heterogeneity.all")){
    dir.create("./00.heterogeneity.all")
  }
  
  
}


library(clusterProfiler)
# Candidated_genes<-read.csv('../../06_lasso_svm.rfe/03.Intersection_gene4.csv')
Candidated_genes <- read.csv('../../02_Hubgene/02.hubgenes.csv',row.names = 1)
id<-bitr(Candidated_genes$symbol,fromType = "SYMBOL",
         toType = c("ENSEMBL"),
         OrgDb = 'org.Hs.eg.db')
id$ENSEMBL<-paste('eqtl','a',id$ENSEMBL,sep='-')
ids<-id$ENSEMBL
i <- a[1]

# for (i in 1:length(ids))  {
for (i in a)  {
  tryCatch({
    
    #{i <- 21
    
    setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/00.mr_ivw_mre/")
    dir.create(paste(ids[i],sep='_'))
    setwd(paste(ids[i],sep='_'))
    
    dat <- read.csv(paste0('../../',ids[i],'/03.harmonized_dat_input.csv'))
    
    #result <- generate_odds_ratios(mr_res = mr(dat))  
    result <- generate_odds_ratios(mr_res = mr(dat, method_list = c("mr_ivw_mre")))   
    write.csv(result,file = '05.result_method.csv',row.names=F)
    write.csv(result,file = paste0("../00.result.all/",ids[i],'_',i,'_result_method.csv'),row.names=F)
    
    
    
    library(patchwork)
    library(ggplot2)
    p1 <- mr_scatter_plot(dat,mr_results = mr(dat,method_list = c('mr_ivw')))[[1]] +
      scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.position = 'top',
            legend.direction = 'vertical',
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('01.scatter.pdf',w=12,h=12,family='Times')
    print(p1)
    dev.off()
    
    png('01.scatter.png',w=12,h=12,units='in',res=600,family='Times')
    print(p1)
    dev.off()
    
    #森林图
    res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre"))
    p2<-mr_forest_plot(res_single)[[1]] +
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            # legend.title = element_blank(),
            # legend.position = c(0.1,0.9),
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.title = element_blank(),
            legend.position = 'none',
            legend.direction = 'vertical',
            plot.title=element_text(size=8),
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('02.forest.pdf',w=18,h=20,family='Times')
    print(p2)
    dev.off()
    
    png('02.forest.png',w=18,h=20,units='in',res=600,family='Times')
    print(p2)
    dev.off()
    
    p3<-mr_funnel_plot(res_single)[[1]] +
      scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.position = 'top',
            legend.direction = 'vertical',
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('03.funnel.pdf',w=10,h=11,family='Times')
    print(p3)
    dev.off()
    
    png('03.funnel.png',w=10,h=11,units='in',res=600,family='Times')
    print(p3)
    dev.off()
    
    #异质性检验
    het <- mr_heterogeneity(dat,method_list = c("mr_ivw_mre"))
    het
    write.csv(het,file = '10.heterogeneity.csv',row.names=F)
    write.csv(het,file = paste0("../00.heterogeneity.all/",ids[i],'_',i,'_heterogeneity.csv'),row.names=F)
    
    
    #逐个剔除检验
    single <- mr_leaveoneout(dat, method = mr_ivw_mre)
    p4<-mr_leaveoneout_plot(single)[[1]] +
      theme_classic(base_size = 20)+
      theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.background = element_blank(),
            legend.position = 'none',
            legend.text = element_text(size = 18))+
      theme(legend.background = element_blank(),
            legend.title = element_blank(),
            legend.position = 'none',
            legend.direction = 'vertical',
            plot.title=element_text(size=8),
            legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
    
    pdf('04.leaveoneout.pdf',w=18,h=18,family='Times')
    print(p4)
    dev.off()
    
    png('04.leaveoneout.png',w=18,h=18,units='in',res=600,family='Times')
    print(p4)
    dev.off()
    
    
    
  }, error = function(e){
    print('error')
  },finally = {
    print(ids[i])
    print(i)
  })
}



#手动循环整理表格-----------
rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "04_MR"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)

library(clusterProfiler)
Candidated_genes <- read.csv('../02_Hubgene/02.hubgenes.csv',row.names = 1)

id<-bitr(Candidated_genes$symbol,fromType = "SYMBOL",
         toType = c("ENSEMBL" ),
         OrgDb = 'org.Hs.eg.db')
id$ENSEMBL<-paste('eqtl','a',id$ENSEMBL,sep='-')
ids<-id$ENSEMBL

path <- getwd()
herb <- list.files(path = "./00.result.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.result.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
results_all <- Reduce(rbind,tem)
results_all$id <- results_all$exposure
results_all <- results_all[,c('id','method','nsnp','pval','or','or_lci95','or_uci95')]


path <- getwd()
herb <- list.files(path = "./00.heterogeneity.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.heterogeneity.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
heterogeneity_all <- Reduce(rbind,tem)
heterogeneity_all <- heterogeneity_all[heterogeneity_all$method=='Inverse variance weighted',]
heterogeneity_all$id <- heterogeneity_all$exposure
heterogeneity_all <- heterogeneity_all[,c('id','Q_pval')]
colnames(heterogeneity_all) <- c('id','heterogeneity_pval')

results_all1 <- merge(results_all,heterogeneity_all,by='id')


# path <- getwd()
# herb <- list.files(path = "./00.mr_ivw_mre/00.result.all/")
# tem <- lapply(herb,function(x){
#   tw <- read.csv(paste0("./00.mr_ivw_mre/00.result.all/",x))
#   #tw$herb = gsub('.csv','',x)
#   tw
# })
# results_all <- Reduce(rbind,tem)
# results_all$id <- results_all$exposure
# results_all <- results_all[,c('id','method','nsnp','pval','or','or_lci95','or_uci95')]
# 
# 
# 
# path <- getwd()
# herb <- list.files(path = "./00.mr_ivw_mre/00.heterogeneity.all/")
# tem <- lapply(herb,function(x){
#   tw <- read.csv(paste0("./00.mr_ivw_mre/00.heterogeneity.all/",x))
#   #tw$herb = gsub('.csv','',x)
#   tw
# })
# heterogeneity_all <- Reduce(rbind,tem)
# heterogeneity_all$id <- heterogeneity_all$exposure
# heterogeneity_all <- heterogeneity_all[,c('id','Q_pval')]
# colnames(heterogeneity_all) <- c('id','heterogeneity_pval')
# 
# results_all2 <- merge(results_all,heterogeneity_all,by='id')
# 
# 
# # i <- 2
# for (i in 1:nrow(results_all2)) {
#   a_row_index <- which(results_all1$id == results_all2$id[i] & results_all1$method == "Inverse variance weighted")
#   results_all1[a_row_index, ] <- results_all2[i, ]
# }


setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
path <- getwd()
herb <- list.files(path = "./00.pleiotropy.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.pleiotropy.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
pleiotropy.all <- Reduce(rbind,tem)
pleiotropy.all$id <- pleiotropy.all$exposure
pleiotropy.all <- pleiotropy.all[,c('id','pval')]
colnames(pleiotropy.all) <- c('id','pleiotropy_pval')

setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
path <- getwd()
herb <- list.files(path = "./00.steiger_dat.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.steiger_dat.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
steiger_dat.all <- Reduce(rbind,tem)
steiger_dat.all$id <- steiger_dat.all$exposure
#steiger_dat.all <- steiger_dat.all[,c('id','SNP','steiger_dir','steiger_pval')]
steiger_dat.all <- steiger_dat.all[,c('id','correct_causal_direction','steiger_pval')]




setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
path <- getwd()
herb <- list.files(path = "./00.R&F.all/")
tem <- lapply(herb,function(x){
  tw <- read.csv(paste0("./00.R&F.all/",x))
  #tw$herb = gsub('.csv','',x)
  tw
})
RF_all <- Reduce(rbind,tem)
RF_all$id <- RF_all$exposure
RF_all <- RF_all[,c('id','SNP','mr_keep','F')]

result.write1 <- results_all1
result.write1 <- merge(result.write1, pleiotropy.all, by='id')
result.write1 <- merge(result.write1, steiger_dat.all, by='id')
result.write1 <- merge(result.write1, RF_all, by='id')

result.write1$ENSEMBL<-substr(result.write1$id,1,22)
result.write1<-merge(result.write1,id,by='ENSEMBL')
write.csv(result.write1,'01.results_all.csv',row.names = F)

length(result.write1$SYMBOL %>% unique())
# 15



#基因挑选--------------
getwd()
rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/04_MR/")
# if (! dir.exists("./09_MR_ieu32_kb10")){
#   dir.create("./09_MR_ieu32_kb10")
# }
# setwd("./09_MR_ieu32_kb10")

# cox <- read.csv('sig_diff_write_25_GSE16879.csv')
# cox <- cox[,c(1,2)]

mr <- read.csv('01.results_all.csv')
table(mr$mr_keep)
mr[order(mr$F),]$F[1]


#查看5种算法
a <- mr[,c(1:6,16),] %>% unique()
a <- a[a$pval<0.05,]

#查看heterogeneity
a <- mr[(mr$method=="Inverse variance weighted" | mr$method=="Inverse variance weighted (multiplicative random effects)") & mr$pval<0.05,]
a <- a[,c(2,3,9,16)] %>% unique()
a <- a[a$heterogeneity_pval<0.05,]

#查看steiger
a <- mr[(mr$method=="Inverse variance weighted" | mr$method=="Inverse variance weighted (multiplicative random effects)") & mr$pval<0.05,]
table(a$correct_causal_direction)
a <- a[a$correct_causal_direction=='TRUE',]
a <- a[as.numeric(a$steiger_pval)<0.05,]

#查看pleiotropy
a <- mr[(mr$method=="Inverse variance weighted" | mr$method=="Inverse variance weighted (multiplicative random effects)") & mr$pval<0.05,]
a <- a[as.numeric(a$pleiotropy_pval)>0.05,]
gene <- a[,c(2,3,16,5,6)] %>% unique()
nrow(gene)
# 10
write.csv(gene,'02.MR_gene.csv',row.names = F)





# 
# rm(list = ls())
# setwd("/data/nas1/liky/project/35_LZZK-10415-3/")
# if (! dir.exists("./03_MR.kb100")){
#   dir.create("./03_MR.kb100")
# }
# setwd("./03_MR.kb100")
# result.write2 <- read.csv('./01.results0.05.csv',row.names = 1)
# gene <- read.csv('./03.uniCox_MR_gene.csv',row.names = 1)
# result.write2 <- result.write2[result.write2$SYMBOL %in% gene$id,]
# out_multi <- unique(result.write2[,c(1:8,16)])
# out_multi <- out_multi[order(out_multi$SYMBOL),]
# 
# #write.csv(out_multi$SYMBOL,'02.MR_gene.csv')
# # out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
# # out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
# hz <- paste(round(out_multi$or,5),
#             "(",round(out_multi$or_lci95,5),
#             "-",round(out_multi$or_uci95,5),")",sep = "")
# 
# 
# tabletext <- cbind(c(NA,"Gene",out_multi$SYMBOL),
#                    c(NA,"P value",round(out_multi$pval,5)),
#                    c(NA,"OR(95% CI)",hz))
# 
# library(forestplot)
# pforest <- forestplot(labeltext=tabletext,
#                       graph.pos=2,  #为Pvalue箱线图所在的位置
#                       col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
#                       mean=c(NA,NA,out_multi$or),
#                       lower=c(NA,NA,out_multi$or_lci95), #95%置信区间下限
#                       upper=c(NA,NA,out_multi$or_uci95), #95%置信区间上限
#                       # xlab="Unicox−analytic Hazard Ratio of immune cell population",
#                       boxsize=0.1,lwd.ci=2,   #箱子大小，线的宽度
#                       ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
#                       zero=1,lwd.zero=1,      #zero线宽 基准线的位置
#                       colgap=unit(5,"mm"),    #列间隙
#                       xticks = c(1, 1.002), #横坐标刻度
#                       lwd.xaxis=1,            #X轴线宽
#                       lineheight = unit(1,"cm"), #固定行高
#                       graphwidth = unit(.3,"npc"), #图在表中的宽度比例
#                       cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
#                       hrzl_lines=list("2" = gpar(lwd=2, col="black"),
#                                       "3" = gpar(lwd=2, col="black") #第三行顶部加黑线，引号内数字标记行位置
#                                       #"15" = gpar(lwd=2, col="black")
#                       ),#最后一行底部加黑线,"29"中数字为nrow(tabletext)+1
#                       mar=unit(rep(0.1, times = 4), "mm"),#图形页边距
#                       # #fpTxtGp函数中的cex参数设置各个组件的大小
#                       txt_gp=fpTxtGp(label=gpar(cex=1),
#                                      ticks=gpar(cex=1.2),
#                                      xlab=gpar(cex = 1.25),
#                                      title=gpar(cex = 1.2))
# )
# print(pforest)
# pdf('03.forest.MR.pdf',w=8,h=5,family='Times')
# print(pforest)
# dev.off()
# 
# png('03.forest.MR.png',w=8,h=5,units='in',res=600,family='Times')
# print(pforest)
# dev.off()
# 
# 
# rm(list = ls())
# setwd("/data/nas1/liky/project/35_LZZK-10415-3/")
# if (! dir.exists("./03_MR.kb100")){
#   dir.create("./03_MR.kb100")
# }
# setwd("./03_MR.kb100")
# cox <- read.csv('../06_uniCox/01.uniCox_0.01.csv')
# gene <- read.csv('./03.uniCox_MR_gene.csv',row.names = 1)
# result.write2 <- cox[cox$id %in% gene$id,]
# out_multi <- result.write2
# out_multi <- out_multi[order(out_multi$id),]
# 
# #write.csv(out_multi$SYMBOL,'02.MR_gene.csv')
# # out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
# # out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
# hz <- paste(round(out_multi$HR,5),
#             "(",round(out_multi$HR.95L,5),
#             "-",round(out_multi$HR.95H,5),")",sep = "")
# 
# 
# tabletext <- cbind(c(NA,"Gene",out_multi$id),
#                    c(NA,"P value",round(out_multi$pval,5)),
#                    c(NA,"HR(95% CI)",hz))
# 
# library(forestplot)
# pforest <- forestplot(labeltext=tabletext,
#                       graph.pos=2,  #为Pvalue箱线图所在的位置
#                       col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
#                       mean=c(NA,NA,out_multi$HR),
#                       lower=c(NA,NA,out_multi$HR.95L), #95%置信区间下限
#                       upper=c(NA,NA,out_multi$HR.95H), #95%置信区间上限
#                       # xlab="Unicox−analytic Hazard Ratio of immune cell population",
#                       boxsize=0.1,lwd.ci=2,   #箱子大小，线的宽度
#                       ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
#                       zero=1,lwd.zero=1,      #zero线宽 基准线的位置
#                       colgap=unit(5,"mm"),    #列间隙
#                       xticks = c(1, 2.8), #横坐标刻度
#                       lwd.xaxis=1,            #X轴线宽
#                       lineheight = unit(1,"cm"), #固定行高
#                       graphwidth = unit(.3,"npc"), #图在表中的宽度比例
#                       cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
#                       hrzl_lines=list("2" = gpar(lwd=2, col="black"),
#                                       "3" = gpar(lwd=2, col="black") #第三行顶部加黑线，引号内数字标记行位置
#                                       #"15" = gpar(lwd=2, col="black")
#                       ),#最后一行底部加黑线,"29"中数字为nrow(tabletext)+1
#                       mar=unit(rep(0.1, times = 4), "mm"),#图形页边距
#                       # #fpTxtGp函数中的cex参数设置各个组件的大小
#                       txt_gp=fpTxtGp(label=gpar(cex=1),
#                                      ticks=gpar(cex=1.2),
#                                      xlab=gpar(cex = 1.25),
#                                      title=gpar(cex = 1.2))
# )
# print(pforest)
# pdf('03.forest.cox.pdf',w=8,h=5,family='Times')
# print(pforest)
# dev.off()
# 
# png('03.forest.cox.png',w=8,h=5,units='in',res=600,family='Times')
# print(pforest)
# dev.off()
# 
# 
# 
# 
# 
# # 
# # out_multi <- unique(result.write2[,c(1:8,16)])
# # out_multi <- out_multi[order(out_multi$or),]
# # write.csv(out_multi$SYMBOL,'02.MR_gene.csv')
# # # out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
# # hz <- paste(round(out_multi$or,5),
# #             "(",round(out_multi$or_lci95,5),
# #             "-",round(out_multi$or_uci95,5),")",sep = "")
# # 
# # 
# # tabletext <- cbind(c(NA,"Gene",out_multi$SYMBOL),
# #                    c(NA,"P value",round(out_multi$pval,5)),
# #                    c(NA,"Hazard Ratio(95% CI)",hz))
# # 
# # library(forestplot)
# # pforest <- forestplot(labeltext=tabletext,
# #                       graph.pos=2,  #为Pvalue箱线图所在的位置
# #                       col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
# #                       mean=c(NA,NA,out_multi$or),
# #                       lower=c(NA,NA,out_multi$or_lci95), #95%置信区间下限
# #                       upper=c(NA,NA,out_multi$or_uci95), #95%置信区间上限
# #                       # xlab="Unicox−analytic Hazard Ratio of immune cell population",
# #                       boxsize=0.3,lwd.ci=2,   #箱子大小，线的宽度
# #                       ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
# #                       zero=1,lwd.zero=1,      #zero线宽 基准线的位置
# #                       colgap=unit(5,"mm"),    #列间隙
# #                       xticks = c(0.998, 1.002), #横坐标刻度
# #                       lwd.xaxis=1,            #X轴线宽
# #                       lineheight = unit(1,"cm"), #固定行高
# #                       graphwidth = unit(.3,"npc"), #图在表中的宽度比例
# #                       cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
# #                       hrzl_lines=list("2" = gpar(lwd=2, col="black"),
# #                                       "3" = gpar(lwd=2, col="black") #第三行顶部加黑线，引号内数字标记行位置
# #                                       #"15" = gpar(lwd=2, col="black")
# #                       ),#最后一行底部加黑线,"29"中数字为nrow(tabletext)+1
# #                       mar=unit(rep(0.1, times = 4), "mm"),#图形页边距
# #                       # #fpTxtGp函数中的cex参数设置各个组件的大小
# #                       txt_gp=fpTxtGp(label=gpar(cex=1),
# #                                      ticks=gpar(cex=1.2),
# #                                      xlab=gpar(cex = 1.25),
# #                                      title=gpar(cex = 1.2))
# # )
# # print(pforest)
# # pdf('03.forest.pdf',w=10,h=8,family='Times')
# # print(pforest)
# # dev.off()
# # 
# # png('03.forest.png',w=10,h=8,units='in',res=600,family='Times')
# # print(pforest)
# # dev.off()
# # 
# 
# 
# 
# # #多因素--------
# # rm(list = ls())
# # setwd("/data/nas1/liky/project/35_LZZK-10415-3/")
# # if (! dir.exists("./03_MR.kb100")){
# #   dir.create("./03_MR.kb100")
# # }
# # setwd("./03_MR.kb100")
# # 
# # if (!dir.exists("02.MVMR")) {dir.create("02.MVMR")}
# # setwd("02.MVMR")
# # 
# # id_exposure <- read.csv('../01.results0.05.csv')
# # id_exposure <- id_exposure$ENSEMBL
# # id_exposure <- unique(id_exposure)
# # 
# # # #在线
# # # library(TwoSampleMR)
# # # #id_exposure <- c('eqtl-a-ENSG00000102760','eqtl-a-ENSG00000116161','eqtl-a-ENSG00000163599')
# # # id_outcome <- "ebi-a-GCST003484"
# # # exposure_dat <- mv_extract_exposures(id_exposure,clump_kb = 10,clump_r2 = 0.001,pval_threshold = 5e-08)
# # # write.csv(exposure_dat,'01.exposure_dat_input.csv',quote=F,row.names=F)
# # # outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) #提取结局数据
# # # write.csv(outcome_dat,'02.outcome_dat_input.csv',quote=F,row.names=F)
# # id_exposure <- ieugwasr::legacy_ids(id_exposure)
# # #本地
# # ##暴露因素
# # exposure_dat <- data.frame()
# # for (i in 1:length(id_exposure))  {
# #   i <- 1
# #   exposure <- extract_instruments(outcomes = id_exposure[i],p1 = 5e-08,r2 = 0.001,kb=10)
# #   vcf1 <- VariantAnnotation::readVcf(paste('/data/nas1/luchunlin/epQTL_coloc/eQTL_vcf/',id_exposure[i],'.vcf.gz',sep=''),"hg19")
# #   # 将vcf转化为TwoSampleMR的格式
# #   exp_data <- gwasvcf_to_TwoSampleMR(vcf1)
# #   
# #   temp_dat <- exp_data[exp_data$pval.exposure<5e-08,] #暴露因素和SNP强相关
# #   #temp_dat <- exp_data
# #   # 补充ld_clump需要的三列
# #   temp_dat$id <- temp_dat$id.exposure # id列用于区分是单个或多个来源的gwas数据，clump按照单个gwas数据进行
# #   temp_dat$rsid <- temp_dat$SNP
# #   temp_dat$pval <- temp_dat$pval.exposure
# #   # 执行ld_clump
# #   # bfile指定参考文件的路径
# #   exp_dat <- ld_clump(temp_dat,
# #                       plink_bin = get_plink_exe(),
# #                       bfile = '/data/nas1/tancj/database/MR/g1000_eur/g1000_eur' ,
# #                       clump_kb = 10, clump_r2 = 0.001)
# #   exp_dat$id.exposure<-id_exposure[i]
# #   exp_dat$exposure = paste(id_exposure[i],' || id:',id_exposure[i],sep='')
# #   #exp_dat <- exp_dat[exp_dat$mr_keep.exposure=='TRUE',]
# #   exposure_dat <- rbind(exposure_dat,exp_dat)
# # }
# # save.image('exposure_dat.RData')
# # load('exposure_dat.RData')
# # 
# # a <- table(exposure_dat$SNP) %>% as.data.frame()
# # a <- a[a$Freq==13,]
# # 
# # b <- exposure_dat[exposure_dat$pval.exposure<5e-08,]
# # b <- b$SNP %>% unique()
# # 
# # exposure_dat2 <- exposure_dat
# # exposure_dat2 <- exposure_dat2$SNP 
# # exposure_dat2 <- exposure_dat2 %>% unique()
# # 
# # exposure_dat2 <- subset(exposure_dat2, exposure_dat2$SNP%in%a$Var1)
# # exposure_dat2 <- subset(exposure_dat2, exposure_dat2$SNP%in%b)
# # 
# # exp_dat <- exp_dat[,c('SNP','exposure','id.exposure','effect_allele.exposure','other_allele.exposure',
# #                       'eaf.exposure','beta.exposure','se.exposure','pval.exposure')]
# # 
# # write.csv(exposure_dat,'01.exposure_dat_input.csv',quote=F,row.names=F)
# # 
# # #结局
# # load('/data/nas1/liky/project/35_LZZK-10415-3/03_MR.kb100/outcome.RData')
# # outcome <- gwasvcf_to_TwoSampleMR(outcome)
# # 
# # #VCF文件
# # outcome_dat <- format_data(
# #   dat=outcome,
# #   type = "outcome",
# #   snps = exposure_dat$SNP,
# #   header = TRUE,
# #   snp_col = "SNP",
# #   beta_col = "beta.exposure",
# #   se_col = "se.exposure",
# #   effect_allele_col = "effect_allele.exposure",
# #   other_allele_col = "other_allele.exposure",
# #   pval_col = "pval.exposure",
# #   ncase_col = "ncase.exposure",
# #   ncontrol_col = "ncontrol.exposure",
# #   chr_col = "chr.exposure",
# #   pos_col = "pos.exposure"
# # )
# # ## ## 解决变量 vcf文件  
# # outcome_dat<- subset(outcome, outcome$SNP%in%exposure_dat$SNP)
# # outcome_dat$id.outcome <- "UC"
# # outcome_dat$outcome <- "UC || id:ebi-a-GCST003484"
# # write.csv(outcome_dat,'02.outcome_dat_input.csv',quote=F,row.names=F)
# # 
# # 
# # 
# # mvdat <- mv_harmonise_data(exposure_dat,  outcome_dat) # 对数据进行harmonisation
# # 
# # 
# # mvdharmonisationres <- mv_multiple(mvdat)
# # mv.result <- mvdharmonisationres$result%>%as.data.frame()
# # mv.result <- generate_odds_ratios(mv.result)
# # mv.result$pval <- round(mv.result$pval,digits = 5)
# # # write.table(mv.result,file = 'MVMR.result.xls',sep = '\t',row.names = F,quote = F)
# # write.table(mv.result,file = '03.MVMR.result.xls',sep = '\t',row.names = F,quote = F)
# # res <- mv_multiple(mvdat,plots=T)
# # 
# # ##画一下森林图------
# # 
# # or <- paste(round(mv.result$or,5),
# #             "(",round(mv.result$or_lci95,5),
# #             "-",round(mv.result$or_uci95,5),")",sep = "")
# # tabletext <- cbind(c(NA,'Exposure',c('RGCC','CACYBP','CTLA4')),
# #                    c(NA,"P value",ifelse(mv.result$pval<0.001,
# #                                          "< 0.001",
# #                                          round(mv.result$pval,5))),
# #                    c(NA,"Odd Ratio(95% CI)",or))
# # 
# # 
# # library(forestplot)
# # pdf(file = "04.mvmr.forest.pdf", height = 6, width =10, onefile = F,family='Times')
# # forestplot(labeltext=tabletext,
# #            graph.pos=4,  #为Pvalue箱线图所在的位置
# #            is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
# #            col=fpColors(box="red", lines="royalblue", zero = "gray50"),
# #            mean=c(NA,NA,mv.result$or),
# #            lower=c(NA,NA,mv.result$or_lci95), #95%置信区间下限
# #            upper=c(NA,NA,mv.result$or_uci95), #95%置信区间上限
# #            boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
# #            ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
# #            zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
# #            colgap=unit(5,"mm"),    #列间隙
# #            xticks = c(0.98,1.02), #横坐标刻度
# #            lwd.xaxis=2,            #X轴线宽
# #            lineheight = unit(1.2,"cm"), #固定行高
# #            graphwidth = unit(.5,"npc"), #图在表中的宽度比例
# #            cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
# #            hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
# #            # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
# #            #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
# #            #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
# #            # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
# #            #fpTxtGp函数中的cex参数设置各个组件的大小
# #            txt_gp=fpTxtGp(label=gpar(cex=1),
# #                           ticks=gpar(cex=0.8, fontface = "bold"),
# #                           xlab=gpar(cex = 1, fontface = "bold"),
# #                           title=gpar(cex = 1.25, fontface = "bold")),
# #            xlab="Odd Ratio",
# #            grid = T) # 垂直于x轴的网格线，对应每个刻度
# # dev.off()
# # 
# # png(file = "04.mvmr.forest.png", height = 6, width = 10, units='in',res=600,family='Times')
# # forestplot(labeltext=tabletext,
# #            graph.pos=4,  #为Pvalue箱线图所在的位置
# #            is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
# #            col=fpColors(box="red", lines="royalblue", zero = "gray50"),
# #            mean=c(NA,NA,mv.result$or),
# #            lower=c(NA,NA,mv.result$or_lci95), #95%置信区间下限
# #            upper=c(NA,NA,mv.result$or_uci95), #95%置信区间上限
# #            boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
# #            ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
# #            zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
# #            colgap=unit(5,"mm"),    #列间隙
# #            xticks = c(0.98,1.02), #横坐标刻度
# #            lwd.xaxis=2,            #X轴线宽
# #            lineheight = unit(1.2,"cm"), #固定行高
# #            graphwidth = unit(.5,"npc"), #图在表中的宽度比例
# #            cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
# #            hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
# #            # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
# #            #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
# #            #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
# #            # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
# #            #fpTxtGp函数中的cex参数设置各个组件的大小
# #            txt_gp=fpTxtGp(label=gpar(cex=1),
# #                           ticks=gpar(cex=0.8, fontface = "bold"),
# #                           xlab=gpar(cex = 1, fontface = "bold"),
# #                           title=gpar(cex = 1.25, fontface = "bold")),
# #            xlab="Odd Ratio",
# #            grid = T) # 垂直于x轴的网格线，对应每个刻度
# # dev.off()




#YT------------------------------------------------------------------
# rm(list = ls())
# setwd("/data/nas1/liky/project/38_BJTC-574-2/")
# if (! dir.exists("./07_MR_kb10/kb_1000")){
#   dir.create("./07_MR_kb10/kb_1000")
# }
# setwd("./07_MR_kb10/kb_1000")
# 
# library(clusterProfiler)
# # library(ClumpMR)
# library(plinkbinr)
# library(VariantAnnotation)
# library(gwasglue)
# library(data.table )
# library(tidyr )
# library(dplyr)
# library(TwoSampleMR)
# library(ieugwasr)
# 
# 
# 
# ###基因
# Candidated_genes<-read.csv('jj_gene25_GSE16879.csv')
# 
# id<-bitr(Candidated_genes$gene,fromType = "SYMBOL",
#          toType = c("ENSEMBL" ),
#          OrgDb = 'org.Hs.eg.db')
# id$ENSEMBL<-paste('eqtl','a',id$ENSEMBL,sep='-')
# ids<-id$ENSEMBL
# 
# 
# exp_data <- readRDS("/data/nas1/tancj/MR_project/01_XA-0224-5/exp_data_2023.4.3_p5e8_idALL.rds")
# 
# 
# ids <- ids[ids %in% unique(exp_data$id.exposure)]
# vcf2 <- VariantAnnotation::readVcf('ieu-a-970.vcf.gz', "hg19")
# outcome <- gwasvcf_to_TwoSampleMR(vcf2)
# outcome_all <- outcome
# # load("../00_Rawdata/outcome_all.Rdata")
# 
# 
# dat.res <- list()
# expo.res <- list()
# results_all<-data.frame()
# pelpvalue <- data.frame()
# options(timeout = 100000)    #等待时间
# i =1
# for (i in 1:length(ids))  {
#   ####暴露因素
#   exp_dat <- exp_data[exp_data$id.exposure == ids[i],]
#   ## 结局变量
#   outcome_dat<- subset(outcome_all,outcome_all$SNP %in% exp_dat$SNP)
#   if (nrow(outcome_dat) < 3) {
#     next
#   }
#   dat <- harmonise_data(exposure_dat =exp_dat,outcome_dat = outcome_dat)
#   ###判断F值---
#   clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
#   clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
#   dat = clumping[clumping$F >10,]
#   ###判断方向---
#   # steiger_sl<-directionality_test(dat)
#   # steiger_sl<-steiger_filtering(dat)
#   # steiger_sl
#   dat.res[[i]] <- dat
#   result <- generate_odds_ratios(mr_res = mr(dat))
#   if (result$nsnp[1] < 3) {
#     next
#   }
#   results_all<-rbind(results_all,result)
#   ple <- mr_pleiotropy_test(dat)
#   pelpvalue <- rbind(pelpvalue,ple)
#   # clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
#   # clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
#   # clumping
#   i = i+1
# }
# 
# # names(dat.res) <- ids
# system.time(save(dat.res, file = "dat.res.Rdata"))
# system.time(save(expo.res, file = "expo.res.Rdata"))
# 
# save.image("MR.all.RData")
# 
# write.csv(results_all,'01.results_all_TrMGs.csv')
# results<-results_all[results_all$method=='Inverse variance weighted',]
# results<-results[results$pval<0.05,]
# results$ENSEMBL<-results$id.exposure
# results<-merge(results,id,by='ENSEMBL')
# write.csv(results,'02.results.csv')
# 
# pelpvalue <- pelpvalue[pelpvalue$pval>0.05,]
# pelpvalue$ENSEMBL <- pelpvalue$id.exposure
# pelpvalue <- merge(pelpvalue,id,by='ENSEMBL')
# 
# gene <- intersect(results$SYMBOL,pelpvalue$SYMBOL)
# write.csv(gene,file="gene.csv")
# 
# 
# 
# 
# # ###-----------------
# # library(ggsci)
# # # devtools::install_github("rondolab/MR-PRESSO")
# # library(MRPRESSO)
# # load("dat.res.Rdata")
# # gene <- read.csv('../00_DrugTarget/tcmsp_herb_gene.csv')
# # id<-bitr(gene$x,fromType = "SYMBOL",
# #          toType = c("ENSEMBL" ),
# #          OrgDb = 'org.Hs.eg.db')
# # id$ENSEMBL<-paste('eqtl','a',id$ENSEMBL,sep='-')
# # ids<-id$ENSEMBL
# # names(dat.res) <- ids[1:205]
# # 
# # 
# # gene <- read.csv("02.results.csv")
# # gene<-merge(gene,id,by='ENSEMBL')
# # dat_all_gene <- dat.res[gene$id.exposure]
# # 
# # 
# # 
# # # #异质性检验+多效性检验+方向性检验------------------
# # heterogeneity_result<-data.frame()
# # pleiotropy_result <- data.frame()
# # out_result <- data.frame()
# # for (i in c(1:length(dat_all_gene))) {
# #   dat <- dat_all_gene[[i]]
# #   sy <- gene$SYMBOL[gene$id.exposure == dat$id.exposure[1]]
# #   het <- mr_heterogeneity(dat)
# #   het$symbol <- sy
# #   pre <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
# #                    OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
# #                    SignifThreshold = 0.05)
# #   pre1 <-pre[[2]]$`Global Test`$Pvalue
# #   ple <- mr_pleiotropy_test(dat)
# #   ple$symble <- sy  
# #   ple$mr_presso <- pre1
# #   out <- directionality_test(dat)
# #   out$symble <- sy
# #   heterogeneity_result<-rbind(heterogeneity_result,het)
# #   pleiotropy_result <- rbind(pleiotropy_result,ple)
# #   out_result <- rbind(out_result,out)
# #   
# # }
# # 
# # write.csv(heterogeneity_result,file = 'heterogeneity.csv',row.names=F)
# # write.csv(pleiotropy_result,file = 'pleiotropy.csv',row.names=F)
# # write.csv(out_result,file = 'Steiger.csv',row.names=F)
# # 
# # 



