rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "14_IC50"
if( !dir.exists(dir)){dir.create(dir)}
setwd(dir)
 
###---- load packages ----####
library(pRRophetic); library(ggplot2); library(rstatix); library(tidyverse); library(ggpubr)
library(ggsci); library(ggplot2); library(ggpol); library(Rmisc)

###---- load packages ----####
# read train expression data file
data <- read.csv(file = "/data/nas3/yuanshilin/project/09.YQFZ-10802-3/00_rawdata/SKCM/01.dat.fpkm.csv", header=TRUE, row.names=1, check.names = FALSE)

# read risk data file
raw.risk <- read.csv("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/07_Riskmodel/01.train_riskScore.csv", row.names = 1) %>% dplyr::select(c("id", "risk"))

# read drug data file 
drug <- read.table(file = '/data/nas1/yuanshilin/Pipline/Medicinal_Sensity/drugs.txt', sep = '\t', header = F)

# check sample 
data <- data[, raw.risk$id]
identical(colnames(data), raw.risk$id)

# assignment
ic.dat <- data

###---- calculated IC50 ----####
# create a empty table
empty <- data.frame(row.names = colnames(data), colnames(data)) 

cnt<-1
while (cnt < 139) {
  print(cnt)
  predictedPtype <- pRRopheticPredict(as.matrix(ic.dat), drug[cnt, ], selection = 1)
  Tipifarnib <- data.frame(predictedPtype)
  colnames(Tipifarnib) <- drug[cnt, ]
  empty <- cbind(empty, Tipifarnib)
  cnt <- cnt + 1
}

# save results
save.image('IC50.RData')
write.csv(empty, '01.IC50.csv')


###---- wilcoxon test ----####
load('IC50.RData')

# read ic50 results 
ic50.res <- read.csv("01.IC50.csv", check.names = F, row.names = 1)

# organize risk data 
risk <- data.frame(id = raw.risk$id, risk = raw.risk$risk)

# check sample
ic50.res <- ic50.res[risk$id, ]

# merge wilcoxon data
wilcox.dat <- merge(ic50.res, risk, by.x = "colnames.data.", by.y = 'id')

# long to short
wilcox.dat <- gather(wilcox.dat, key = drug, value = IC50, -c('colnames.data.', 'risk')) %>% data.frame()

# wilcoxon test 
stat.test <- wilcox.dat %>%
  group_by(drug) %>%
  wilcox_test(IC50 ~ risk) %>%
  adjust_pvalue(method = 'fdr')

# add significant sign
stat.test$p1 <- ifelse(stat.test$p < 0.0001,"****",
                       ifelse(stat.test$p < 0.001,"***",
                              ifelse(stat.test$p < 0.01,"**",
                                     ifelse(stat.test$p < 0.05,"*",'ns'))))

# extract result
ic50.diff <- ic50.res[,stat.test$drug]
identical(stat.test$drug, colnames(ic50.diff))
stat.test$direction <- NA

# delete zero
ic50.diff[is.na(ic50.diff)] <- 0

# add group
for (i in 1:138) {
  stat.test$direction[i] <- ifelse(stat.test$p[i] > 0.05,'no sense',
                                   ifelse(median(as.numeric(ic50.diff[risk$risk=='High',][,i])) > median(as.numeric(ic50.diff[risk$risk=='Low',][,i])),
                                          'Sensitive to low-risk group','Sensitive to high-risk group'))
}

# save wilcoxon test all result 
write.csv(stat.test,'02.IC50_all_wilcox.csv')

# save different drug ic50
stat.test.diff <- stat.test[stat.test$p < 0.05,]
write.csv(stat.test.diff, '03.IC50_diff_wilcox.csv')
paste0("共获得",nrow(stat.test.diff),"个在组间具有差异的化疗药物。")
# "共获得81个在组间具有差异的化疗药物。"


###---- ic50 fig data ----####
table(stat.test$direction)
# no sense Sensitive to high-risk group  Sensitive to low-risk group 
# 57                           39                           42

# organize red drug
red <- data.frame(drug = stat.test$drug[stat.test$direction == "Sensitive to high-risk group"])
# add x and y
red$x <- c(rep(1, 21), rep(2, 18))
red$y <- c(seq(1, 21, 1), seq(4, 21, 1))

# save red data file
write.csv(red, "04.plot_red.csv", row.names = F)

# organize green drug
green <- data.frame(drug = stat.test$drug[stat.test$direction == "Sensitive to low-risk group"])
# add x and y
green$x <- c(rep(1, 21), rep(2, 21))
green$y <- c(seq(1, 21, 1), seq(1, 21, 1))

# save green data file 
write.csv(green, "04.plot_green.csv", row.names = F)

# organize gray drug
gray <- data.frame(drug = stat.test$drug[stat.test$direction == "no sense"])
# add x and y
gray$x <- c(rep(1, 21), rep(2, 21), rep(3, 15))
gray$y <- c(seq(1, 21, 1), seq(1, 21, 1), seq(7, 21, 1))

# save green data file 
write.csv(gray, "04.plot_gray.csv", row.names = F)


###---- IC50 fig ----####
# p0 data 
plot.p0 <- stat.test
num <- table(plot.p0$direction) %>% data.frame()

# organzie p0 data
p0.data <- data.frame(Parties = factor(c("Sensitive to high-risk group", "no sense", "Sensitive to low-risk group"),
                                       levels = c("Sensitive to high-risk group", "no sense", "Sensitive to low-risk group")),
                      seats = c(num$Freq[2], num$Freq[1], num$Freq[3]),   #需要与table的数量一致
                      colors = c("#E56F5E", "grey", "#43978F"),
                      stringsAsFactors = FALSE)

# p0 fig
p0 <- ggplot(p0.data) + 
  geom_parliament(aes(seats = seats, fill = Parties), color = "black") + 
  scale_fill_manual(values = p0.data$colors, labels = p0.data$Parties) +
  coord_fixed() + 
  theme_void()+
  ggtitle("Medicinal Sensity")+
  theme(plot.title = element_text(color = 'black', size = 24, face = "bold", hjust = 0.5),
        text = element_text(family = 'Times'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.text = element_text(color = 'black', size = 12, face = "bold"),
        legend.title = element_text(color = 'black', size = 18, face = "bold"),
        legend.position = "top")
p0

# read red data file
data.p1 <- read.csv('04.plot_red.csv')
table(data.p1$x)
range(data.p1$y)

# p1 fig add red table
p1 <- ggplot(data.p1, aes(x = x, y = y)) +
  geom_point(color = 'white') +
  geom_text(label = data.p1$drug, check_overlap = T,
            vjust = 0.5, size = 4, hjust = 0, nudge_x = 0, alpha = 0.9, col = '#E56F5E')+
  theme_void()+
  xlim(1, max(data.p1$x)+.5) +
  ylim(1, max(data.p1$y)) +
  theme(text = element_text(angle = 0, hjust = 0, family = 'Arial'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "#E56F5E", fill = NA, linewidth = 1))
p1

# read gray data file
data.p2 <- read.csv('04.plot_gray.csv')
table(data.p2$x)
range(data.p2$y)

# p2 fig add gray table
p2 <- ggplot(data.p2, aes(x = x, y = y)) +
  geom_point(color = 'white') +
  geom_text(label = data.p2$drug, check_overlap = T,
            vjust = 0.5, size = 4, hjust = 0, nudge_x = 0, alpha = 0.9, col = 'grey') +
  theme_void() +
  xlim(1, max(data.p2$x)+.5) +
  ylim(1, max(data.p2$y)) +
  theme(text = element_text(angle = 0, hjust = 0, family = 'Arial'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1))
p2

ra# read gray data file
data.p3 <- read.csv('04.plot_green.csv')
table(data.p3$x)
range(data.p3$y)

# p3 fig add green table
p3 <- ggplot(data.p3, aes(x = x, y = y)) +
  geom_point(color = 'white') +
  geom_text(label = data.p3$drug, check_overlap = T,
            vjust = 0.5, size = 4, hjust = 0, nudge_x = 0, alpha = 0.9, col='#43978F')+
  theme_void()+
  xlim(1, max(data.p3$x)+.5) +
  ylim(1, max(data.p3$y)) +
  theme(text = element_text(angle = 0, hjust = 0, family='Arial'),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "#43978F", fill = NA, linewidth = 1))
p3

# arrange fig 
p <- ggarrange(p3, p2, # 要排版的图形
               ncol = 2, nrow = 1, # 行数和列数为2
               widths = c(2, 3)) # 宽度为2：1
p
pp <- ggarrange(p, p1, # 要排版的图形
                ncol = 2, nrow = 1, # 行数和列数为2
                widths = c(5, 2)) # 宽度为2：1
pp
ppp <- ggarrange(p0, pp, # 要排版的图形
                 ncol = 1, nrow = 2, # 行数和列数为2
                 widths = c(1, 1)) # 宽度为2：1
ppp

# save ic50 fig
ggsave(filename = "05.IC50.pdf", ppp, width = 16, height = 9, family = 'Times')
ggsave(filename = "05.IC50.png", ppp, width = 16, height = 9, units = "in", dpi = 600, bg = "white")


