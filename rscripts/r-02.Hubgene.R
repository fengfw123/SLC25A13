rm(list = ls())
setwd("/data/nas3/yuanshilin/project/09.YQFZ-10802-3/")
dir <- "02_Hubgene"
if( !dir.exists(dir)){
  dir.create(dir)
}
setwd(dir)

###---- load packages ----####
pacman::p_load(tidyr, dplyr, VennDiagram, ggvenn,ggplot2)

###---- load rawdata  ----####
# different of tumor and normal
dif.dat <- data.frame(symbol = read.csv("../01_DEGs/02.DEGs_sig(GSE156052).csv",row.names = 1) %>% rownames())

# Mitochondrial metabolism
mi.gene <- readxl::read_xlsx("../00_rawdata/MRGs/Mitochondrial_metabolism.xlsx")

# venn list 
venn.ls <- list("mitochondria gene" = mi.gene$id,
                "different gene" = dif.dat$symbol
               )

# venn plot
venn.plot <- function(data){
  color.set <- c("#49afca", "#d64c39", "#F26800")
  name.color.set <- c("#6a84a0","#a6666b","#6B4945")
  venn.p <- ggvenn(data,
                   fill_color = color.set[1:length(data)],
                   show_percentage = T,
                   stroke_alpha = 0.5,
                   stroke_size = 1,
                   text_size = 5,
                   stroke_color = "white",
                   stroke_linetype = "dashed",
                   set_name_color = name.color.set[1:length(data)],
                   set_name_size = 6,
                   text_color = 'black')
  return(venn.p)
}

# output venn
venn <- venn.plot(venn.ls)
venn
ggsave(filename = "01.hubgene_venn.pdf", venn, width = 8, height = 8, bg = "white")
ggsave(filename = "01.hubgene_venn.png", venn, width = 8, height = 8, bg = "white", dpi = 600, units = "in")

# extract hubgene
hubgenes <- data.frame(symbol = intersect(mi.gene$id, dif.dat$symbol))

# hubgene comma-separated values table
write.csv(hubgenes, file = "02.hubgenes.csv", row.names = T, quote = F)

