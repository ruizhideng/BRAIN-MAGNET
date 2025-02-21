library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(readxl)
library(tidyr)
library(data.table)
library(wordcloud)
library(openxlsx)
library(BSgenome.Hsapiens.UCSC.hg38)
library(phastCons100way.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
library(purrr)
library(PRROC)
options(scipen = 999)

setwd("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Scripts/")

#########################
###Fig. 1B and Fig. 2B###
#########################

NSC_rank<- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>%
  select(c("log2(avg.NSC.RPP+1)", "5_Categories")) %>%
  unique() %>% 
  arrange(`log2(avg.NSC.RPP+1)`) %>%
  mutate(id=row_number())
colnames(NSC_rank) <- c("log2(RPP+1)", "5_Categories", "id")
  
ESC_rank <- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 2) %>%
  select(c("log2(avg.ESC.RPP+1)", "5_Categories")) %>%
  unique() %>% 
  arrange(`log2(avg.ESC.RPP+1)`) %>%
  mutate(id=row_number())
colnames(ESC_rank) <- c("log2(RPP+1)", "5_Categories", "id")

### Fig. 1B (NSC activity Rank)
ggplot(NSC_rank, aes(id, `log2(RPP+1)`)) + 
  geom_point(size=0.5, color = "#CC0000") +
  scale_x_continuous(breaks=seq(0,149000,29800), labels = seq(0,10,by=2))+
  geom_vline(xintercept = c(134100), linetype="dashed") +
  xlab("NSCs enhancers") + ylab("log2(RPP + 1)") +
  theme(axis.title = element_text(size=16), axis.text =  element_text(size=16)) 

### Fig. 2B (ESC activity Rank)
ggplot(ESC_rank, aes(id, `log2(RPP+1)`)) + 
  geom_point(size=0.5, color = "#000099") +
  scale_x_continuous(breaks=seq(0,149000,29800), labels = seq(0,10,by=2))+
  geom_vline(xintercept = c(134100), linetype="dashed") +
  xlab("ESCs enhancers") + ylab("log2(RPP + 1)") +
  theme(axis.title = element_text(size=16), axis.text =  element_text(size=16)) +
  ylim(0,8)

rm(list=ls())
gc()

##################################
###Fig. 1D, Fig. 1E and Fig. 1H###
###Fig. 2C, Fig. 2D and Fig. 2E###
##################################

NSC_color <- c("#FFE6E6", "#FFCCCC","#FF9999","#FF6666","#FF3333","#FF0000","#CC0000")
ESC_color <- c("#CCCCFF","#9999FF","#6666FF","#6666FF","#0000FF","#0000CC","#000099")

### 
### Fig. 1D (NSC target gene expression)
### 
NCREs<- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>%
  .[,c(1,5,7,8)] %>% unique() %>%
  separate(col = NCREs, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end))
colnames(NCREs) <- c("chr", "start", "end", "Cateogory", "GeneType", "GeneName")

NCREs_top10<- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>%
  .[,c(1,4,7,8)] %>%
  subset(`10_Categories` == "Category_10") %>%
  separate(col = NCREs, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end))
NCREs_top10[, "10_Categories"] <- "Top10"
colnames(NCREs_top10) <- c("chr", "start", "end", "Cateogory", "GeneType", "GeneName")

NCREs_top10_OMIM<- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>%
  subset(`10_Categories` == "Category_10") %>%
  filter(!is.na(OMIM_Phenotypes)) %>%
  .[,c(1,4,7,8)] %>%
  separate(col = NCREs, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end)) %>% unique()
NCREs_top10_OMIM[, "10_Categories"] <- "Top10_OMIM"
colnames(NCREs_top10_OMIM) <- c("chr", "start", "end", "Cateogory", "GeneType", "GeneName")

NCREs <- rbind(NCREs, NCREs_top10, NCREs_top10_OMIM) %>% unique()

expression <- read_excel("data/Supplementary Table 8 - Gene expression (RNA-seq).xlsx")

NCREs_gene.expression<- NCREs %>%
  merge(expression[,c(2,8,9)], by = c("GeneName")) %>%
  subset(GeneType %in% c("protein_coding", "protein-coding")) %>%
  mutate(ESC=log2(ESC.YY1_mean+1)) %>% 
  mutate(NSC=log2(NSC.YY1_mean+1)) %>% 
  subset(GeneType %in%  c("protein-coding","protein_coding"))%>% unique() %>%
  dplyr :: select ("GeneName","Cateogory" ,"NSC") %>% 
  unique() %>%
  na.omit() %>%
  gather(key = "Cell_Type", value = "value", -GeneName, -Cateogory) %>%
  unique()

ggplot(NCREs_gene.expression, aes(x=Cateogory, y=value, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  ylab("log2(RPKM+1)")+
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text()) +
  scale_fill_manual(name="NCREs", values= NSC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_gene.expression, NCREs_gene.expression$Cateogory == "Category_4")$value, 
            subset(NCREs_gene.expression, NCREs_gene.expression$Cateogory == "Category_5")$value)

rm(NCREs_top10, NCREs_top10_OMIM, NCREs_gene.expression, expression)

### 
### Fig. 1E (NSC target gene pLI score)
### 
pLI <- read.table("data/pLI.score.bed"); colnames(pLI) <- c("GeneName", "pLI")
NCREs_gene.pLI<- NCREs %>%
  merge(pLI, by = c("GeneName")) %>%
  na.omit() %>%
  dplyr::select(c("GeneName", "pLI", "Cateogory")) %>%
  unique()

ggplot(NCREs_gene.pLI, aes(x=Cateogory, y=pLI, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  ylab("log2(RPKM+1)")+
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text()) +
  scale_fill_manual(name="NCREs", values= NSC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_gene.pLI, NCREs_gene.pLI$Cateogory == "Top10")$pLI, 
            subset(NCREs_gene.pLI, NCREs_gene.pLI$Cateogory == "Top10_OMIM")$pLI)

rm(NCREs_gene.pLI, pLI)

### 
### Fig. 1H (NSC NCREs category, sequence features) 
### 

## ncER percentile
ncER <- read.table("data/Scaffolds_nCER_hg38.bed"); colnames(ncER) <- c("chr", "start", "end", "mean", "median")
NCREs_ncER<- NCREs %>%
  left_join(ncER, by = c("chr", "start", "end")) %>%
  na.omit() %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "median")) %>%
  unique()

ggplot(NCREs_ncER, aes(x=Cateogory, y=median, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("ncER Percentile") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = NSC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_ncER, NCREs_ncER$Cateogory == "Category_5")$median, 
            subset(NCREs_ncER, NCREs_ncER$Cateogory == "Top10")$median)


rm(ncER, NCREs_ncER)

## GC content
NCREs_coordinate <- NCREs %>%
  subset(chr %in% c(paste0("chr", 1:22), "chrX")) %>%
  unite("tmp",chr:start, sep = ":") %>%
  unite("coordinate", tmp:end, sep = "-")

getGC <- function(genome, ranges){
  s <- BSgenome::getSeq(genome, ranges)
  a <- Biostrings::alphabetFrequency(s)
  gc <- apply(a[,c("G","C")], 1, sum)
  all <- apply(a[,c("A", "T", "C", "G")], 1, sum)
  ranges$gc.content <- gc/all
  ranges
}

NCREs_GC.tmp <- data.frame(getGC(BSgenome.Hsapiens.UCSC.hg38, GRanges(unique(NCREs_coordinate$coordinate))))
colnames(NCREs_GC.tmp) <- c("chr", "start", "end", "width", "strand", "GC_content")

NCREs_GC <- NCREs %>% 
  merge(NCREs_GC.tmp, by = c("chr", "start", "end")) %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "GC_content")) %>%
  unique()

ggplot(NCREs_GC, aes(x=Cateogory, y=GC_content, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("GC Score") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = NSC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_GC, NCREs_GC$Cateogory == "Category_2")$GC_content, 
            subset(NCREs_GC, NCREs_GC$Cateogory == "Category_1")$GC_content)

rm(NCREs_GC, NCREs_GC.tmp, getGC)

## phastcons Score
NCREs_phastCons.tmp <- as.data.frame(gscores(phastCons100way.UCSC.hg38, GRanges(NCREs_coordinate$coordinate)))
colnames(NCREs_phastCons.tmp) <- c("chr", "start", "end", "width", "strand", "phastcons")

NCREs_phastCons <- NCREs %>% 
  merge(NCREs_phastCons.tmp, by = c("chr", "start", "end")) %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "phastcons")) %>%
  unique()

ggplot(NCREs_phastCons, aes(x=Cateogory, y=phastcons, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("phastcons Score") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = NSC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_phastCons, NCREs_phastCons$Cateogory == "Category_5")$phastcons, 
            subset(NCREs_phastCons, NCREs_phastCons$Cateogory == "Category_4")$phastcons)

rm(NCREs_phastCons.tmp, NCREs_phastCons, NCREs_coordinate)

## Orion score
orion <- read.table("data/Scaffolds_Orion_hg38.bed"); colnames(orion) <- c("chr", "start", "end", "mean", "median")
NCREs_orion<- NCREs %>%
  left_join(orion, by = c("chr", "start", "end")) %>%
  na.omit() %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "median")) %>%
  unique()

ggplot(NCREs_orion, aes(x=Cateogory, y=median, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("Orion Score") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = NSC_color) +
  ylim(-0.75, 0.2)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_orion, NCREs_orion$Cateogory == "Top10_OMIM")$median, 
            subset(NCREs_orion, NCREs_orion$Cateogory == "Top10")$median)

rm(orion, NCREs_orion)

## CADD PHRED
cadd <- read.table("data/Scaffolds_CADD_hg38.bed"); colnames(cadd) <- c("chr", "start", "end", "mean", "median")
NCREs_cadd <- NCREs %>%
  left_join(cadd, by = c("chr", "start", "end")) %>%
  na.omit() %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "median")) %>%
  unique()

ggplot(NCREs_cadd, aes(x=Cateogory, y=median, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("CADD PHRED") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = NSC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_cadd, NCREs_cadd$Cateogory == "Top10_OMIM")$median, 
            subset(NCREs_cadd, NCREs_cadd$Cateogory == "Top10")$median)

rm(cadd, NCREs_cadd)

### 
### Fig. 2C (ESC target gene expression)
### 
NCREs<- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 2) %>%
  .[,c(1,5,7,8)] %>% unique() %>%
  tidyr::separate(col = NCREs, into = c("chr", "start_tmp"), sep = ":") %>%
  tidyr::separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end))
colnames(NCREs) <- c("chr", "start", "end", "Cateogory", "GeneType", "GeneName")

NCREs_top10<- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 2) %>%
  .[,c(1,4,7,8)] %>%
  subset(`10_Categories` == "Category_10") %>%
  separate(col = NCREs, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end))
NCREs_top10[, "10_Categories"] <- "Top10"
colnames(NCREs_top10) <- c("chr", "start", "end", "Cateogory", "GeneType", "GeneName")

NCREs_top10_OMIM<- read_excel("data/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 2) %>%
  subset(`10_Categories` == "Category_10") %>%
  filter(!is.na(OMIM_Phenotypes)) %>%
  .[,c(1,4,7,8)] %>%
  separate(col = NCREs, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end)) %>% unique()
NCREs_top10_OMIM[, "10_Categories"] <- "Top10_OMIM"
colnames(NCREs_top10_OMIM) <- c("chr", "start", "end", "Cateogory", "GeneType", "GeneName")

NCREs <- rbind(NCREs, NCREs_top10, NCREs_top10_OMIM) %>% unique()

expression <- read_excel("data/Supplementary Table 8 - Gene expression (RNA-seq).xlsx")

NCREs_gene.expression<- NCREs %>%
  merge(expression[,c(2,8,9)], by = c("GeneName")) %>%
  subset(GeneType %in% c("protein_coding", "protein-coding")) %>%
  mutate(ESC=log2(ESC.YY1_mean+1)) %>% 
  mutate(NSC=log2(NSC.YY1_mean+1)) %>% 
  subset(GeneType %in%  c("protein-coding","protein_coding"))%>% unique() %>%
  dplyr :: select ("GeneName","Cateogory" ,"ESC") %>% 
  unique() %>%
  na.omit() %>%
  gather(key = "Cell_Type", value = "value", -GeneName, -Cateogory) %>%
  unique()

ggplot(NCREs_gene.expression, aes(x=Cateogory, y=value, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  ylab("log2(RPKM+1)")+
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text()) +
  scale_fill_manual(name="NCREs", values= ESC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_gene.expression, NCREs_gene.expression$Cateogory == "Top10")$value, 
            subset(NCREs_gene.expression, NCREs_gene.expression$Cateogory == "Top10_OMIM")$value)

rm(NCREs_top10, NCREs_top10_OMIM, NCREs_gene.expression, expression)

### 
### Fig. 2D (ESC target gene pLI score)
### 
pLI <- read.table("data/pLI.score.bed"); colnames(pLI) <- c("GeneName", "pLI")
NCREs_gene.pLI<- NCREs %>%
  merge(pLI, by = c("GeneName")) %>%
  na.omit() %>%
  dplyr::select(c("GeneName", "pLI", "Cateogory")) %>%
  unique()

ggplot(NCREs_gene.pLI, aes(x=Cateogory, y=pLI, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  ylab("log2(RPKM+1)")+
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text()) +
  scale_fill_manual(name="NCREs", values= ESC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_gene.pLI, NCREs_gene.pLI$Cateogory == "Category_5")$pLI, 
            subset(NCREs_gene.pLI, NCREs_gene.pLI$Cateogory == "Top10")$pLI)

rm(NCREs_gene.pLI, pLI)

### 
### Fig. 2E (ESC NCREs category, sequence features) 
### 

## ncER percentile
ncER <- read.table("data/Scaffolds_nCER_hg38.bed"); colnames(ncER) <- c("chr", "start", "end", "mean", "median")
NCREs_ncER<- NCREs %>%
  left_join(ncER, by = c("chr", "start", "end")) %>%
  na.omit() %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "median")) %>%
  unique()

ggplot(NCREs_ncER, aes(x=Cateogory, y=median, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("ncER Percentile") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = ESC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_ncER, NCREs_ncER$Cateogory == "Top10")$median, 
            subset(NCREs_ncER, NCREs_ncER$Cateogory == "Category_5")$median)

rm(ncER, NCREs_ncER)

## GC content
NCREs_coordinate <- NCREs %>%
  subset(chr %in% c(paste0("chr", 1:22), "chrX")) %>%
  unite("tmp",chr:start, sep = ":") %>%
  unite("coordinate", tmp:end, sep = "-")

getGC <- function(genome, ranges){
  s <- BSgenome::getSeq(genome, ranges)
  a <- Biostrings::alphabetFrequency(s)
  gc <- apply(a[,c("G","C")], 1, sum)
  all <- apply(a[,c("A", "T", "C", "G")], 1, sum)
  ranges$gc.content <- gc/all
  ranges
}

NCREs_GC.tmp <- data.frame(getGC(BSgenome.Hsapiens.UCSC.hg38, GRanges(unique(NCREs_coordinate$coordinate))))
colnames(NCREs_GC.tmp) <- c("chr", "start", "end", "width", "strand", "GC_content")

NCREs_GC <- NCREs %>% 
  merge(NCREs_GC.tmp, by = c("chr", "start", "end")) %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "GC_content")) %>%
  unique()

ggplot(NCREs_GC, aes(x=Cateogory, y=GC_content, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("GC Score") +
  ylim(0,1) + 
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = ESC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_GC, NCREs_GC$Cateogory == "Top10")$GC_content, 
            subset(NCREs_GC, NCREs_GC$Cateogory == "Category_5")$GC_content)


rm(NCREs_GC, NCREs_GC.tmp, getGC)

## phastcons Score
NCREs_phastCons.tmp <- as.data.frame(gscores(phastCons100way.UCSC.hg38, GRanges(NCREs_coordinate$coordinate)))
colnames(NCREs_phastCons.tmp) <- c("chr", "start", "end", "width", "strand", "phastcons")

NCREs_phastCons <- NCREs %>% 
  merge(NCREs_phastCons.tmp, by = c("chr", "start", "end")) %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "phastcons")) %>%
  unique()

ggplot(NCREs_phastCons, aes(x=Cateogory, y=phastcons, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("phastcons Score") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = ESC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_phastCons, NCREs_phastCons$Cateogory == "Top10")$phastcons, 
            subset(NCREs_phastCons, NCREs_phastCons$Cateogory == "Top10_OMIM")$phastcons)

rm(NCREs_phastCons.tmp, NCREs_phastCons, NCREs_coordinate)

## Orion score
orion <- read.table("data/Scaffolds_Orion_hg38.bed"); colnames(orion) <- c("chr", "start", "end", "mean", "median")
NCREs_orion<- NCREs %>%
  left_join(orion, by = c("chr", "start", "end")) %>%
  na.omit() %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "median")) %>%
  unique()

ggplot(NCREs_orion, aes(x=Cateogory, y=median, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("Orion Score") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = ESC_color) +
  ylim(-0.75, 0.2)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_orion, NCREs_orion$Cateogory == "Top10_OMIM")$median, 
            subset(NCREs_orion, NCREs_orion$Cateogory == "Top10")$median)

rm(orion, NCREs_orion)

## CADD PHRED
cadd <- read.table("data/Scaffolds_CADD_hg38.bed"); colnames(cadd) <- c("chr", "start", "end", "mean", "median")
NCREs_cadd <- NCREs %>%
  left_join(cadd, by = c("chr", "start", "end")) %>%
  na.omit() %>%
  dplyr::select(c("chr", "start", "end", "Cateogory", "median")) %>%
  unique()

ggplot(NCREs_cadd, aes(x=Cateogory, y=median, fill=Cateogory)) +
  geom_boxplot(lwd=0.35, outlier.colour="black", outlier.shape=16, 
               outlier.size=0.8, alpha=0.8) +
  scale_y_continuous("CADD PHRED") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, color="black", angle = 45, hjust = 1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text())  +
  scale_fill_manual(name="NCREs", values = ESC_color)

# Category_1, Top10, Top10_OMIM
wilcox.test(subset(NCREs_cadd, NCREs_cadd$Cateogory == "Top10_OMIM")$median, 
            subset(NCREs_cadd, NCREs_cadd$Cateogory == "Category_5")$median)


rm(cadd, NCREs_cadd)

#########################
###Fig. 1F and Fig. 2F###
#########################
### Fig. 1F (NSC target gene GO enrichment) --> load sheet_nums <- 1:6 
### Fig. 2D (ESC target gene GO enrichment) --> load sheet_nums <- 7:12

# Define a function to process the data
process_go_terms <- function(file_path, sheet_num, category_name) {
  go_terms <- read.xlsx(file_path, sheet = sheet_num)
  go_terms$Term <- sub("\\s*\\(.*", "", go_terms$Term)
  go_terms$log10pvalue <- (-log10(go_terms$`P-value`))
  go_terms <- go_terms[, c(1, 10)] %>% 
    unique() %>% 
    arrange(desc(log10pvalue)) %>%
    subset(!(Term %in% "regulation of chemokine"))
  go_terms$Term <- gsub("negative regulation", "neg. reg.", go_terms$Term)
  go_terms$Term <- gsub("positive regulation", "pos. reg.", go_terms$Term)
  go_terms$Term <- gsub("regulation", "reg.", go_terms$Term)
  go_terms$Term <- gsub("polymerase", "pol", go_terms$Term)
  
  colnames(go_terms) <- c("Term", category_name)
  return(go_terms)
}

# Process each sheet and create word cloud
create_wordcloud <- function(file_path, sheet_nums, category_names) {
  for (i in seq_along(sheet_nums)) {
    go_terms <- process_go_terms(file_path, sheet_nums[i], category_names[i])
    colors <- brewer.pal(9, "Reds") ### Fig. 1F (NSC target gene GO enrichment)
    # colors <- brewer.pal(9, "Blues") ### Fig. 2D (ESC target gene GO enrichment)
    pdf(paste0("figures/wordcloud_NSC_", category_names[i], ".pdf"))
    wordcloud(words = go_terms$Term, 
              freq = go_terms[[2]], min.freq = 1, max.words = 200,
              random.order = FALSE, rot.per = 0.35, scale = c(2, 0.01),
              colors = colors)
    dev.off()
  }
}

# Define file path, sheet numbers, and category names
file_path <- "data/Supplementary Table 6 - GO BP NSCs and ESCs.xlsx"

sheet_nums <- 1:6 ### Fig. 1F (NSC target gene GO enrichment)
# sheet_nums <- 7:12 ### Fig. 2D (ESC target gene GO enrichment)
category_names <- c("Category_1", "Category_2", "Category_3", "Category_4", "Category_5", "Top10")

# Call the function to create word clouds
create_wordcloud(file_path, sheet_nums, category_names)

rm(list=ls())
gc()


#### TEST: dot plot
file_path <- "data/Supplementary Table 6 - GO BP NSCs and ESCs.xlsx"
go_terms <- read.xlsx(file_path, sheet = 1)
go_terms$Term <- sub("\\s*\\(.*", "", go_terms$Term)
# go_terms$Overlap <- sub("\\s*\\/.*", "", go_terms$Overlap)
go_terms$test <- go_terms$Overlap

go_terms <- go_terms %>% arrange(`Adjusted.P-value`) %>%
  subset(!(Term %in% "regulation of chemokine"))
go_terms$Term <- gsub("negative regulation", "neg. reg.", go_terms$Term)
go_terms$Term <- gsub("positive regulation", "pos. reg.", go_terms$Term)
go_terms$Term <- gsub("regulation", "reg.", go_terms$Term)
go_terms$Term <- gsub("polymerase", "pol", go_terms$Term)



colors <- brewer.pal(9, "Reds") ### Fig. 1F (NSC target gene GO enrichment)

ggplot(go_terms[1:10,]) +
  geom_point(aes(x = Overlap, y = Term, color = `Adjusted.P-value`)) +
  theme_bw() +  
  theme(axis.text.x=element_text(size=9, color="black"),
                      axis.text.y=element_text(size=9, color="black"),
                      axis.title.x=element_blank(),
                      axis.title.y=element_text()) +
  xlab("Gene ratios") +
  ylab("Top 10 significant GO terms") +
  ggtitle("Dotplot of top 10 significant GO terms")


#############
###Fig. 3C###
#############

### Fig. 3C: TF-modisco crucial motifs

excel_file <- "data/Supplementary Table 14 - TF-modisco motifs.xlsx"
all_tabs <- excel_sheets(excel_file)
tabs_to_read <- all_tabs[1:6] 
data_list <- map(tabs_to_read, ~ read_excel(excel_file, sheet = .x) %>%
                   mutate(nor_pval = -log10(pval0)) %>%
                   mutate(tab_name = .x))

# Combine all data frames into a single data frame
combined_data <- bind_rows(data_list, .id = "tab_name")

# Filter data for qval0 < 0.05 and pattern_type == "pos_patterns"
filtered_data <- combined_data %>%
  subset((qval0 < 0.05 & pattern_type == "pos_patterns") | 
           grepl("^SOX", motif_GeneName0) | 
           grepl("^MAZ", motif_GeneName0) |
           grepl("^RFX", motif_GeneName0) |
           grepl("^ERG", motif_GeneName0) |
           grepl("^FLI", motif_GeneName0) |
           grepl("^ETV", motif_GeneName0))

filtered_data <- combined_data %>%
  subset((qval0 < 0.05 & pattern_type == "pos_patterns") )

# Create a custom color palette
my_colors <- c("#234471", "#cb1517")
# Create the plot with facets based on the tab_name
ggplot(combined_data, aes(x = num_seqlets, y = nor_pval)) + 
  geom_point(aes(color = pattern_type, size = num_seqlets), alpha = 0.7) +
  geom_text_repel(data = filtered_data, aes(label = motif_GeneName0), size = 3, max.overlaps = 10) +
  xlab("Number of seqlets") + ylab("-log10(pvalue)")  + 
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ tab_name, scales = "free_x") +
  coord_cartesian(xlim = c(0, 1500)) + # Limit x-axis to 0-1500
  theme(panel.grid.major = element_line(color = "gray", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = "gray", size = 0.5),
        strip.text = element_text(color = "black", size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 8)) 

rm(list=ls())
gc()

#####################
###Fig. 3B and S5C###
#####################

### Fig. 3B and S5C: PRROC for AUPRC in Top10 and Bottom10

preds_targets<- read.csv("data/Top_Bottom_ESC.csv")
# calculate precision and recall values
pr <- pr.curve(scores.class0 = preds_targets$preds, weights.class0 = preds_targets$label, curve = TRUE)
plot(pr)

roc <- roc.curve(scores.class0 = preds_targets$preds, weights.class0 = preds_targets$label, curve = TRUE)
plot(roc)

##############
###Fig. S5B###
##############

### Fig. S5B: preds, targets correlation of BRAIN-MAGNET

library(smplot2)
library(ggrastr)
library(reticulate)
library(MASS)
library(viridis)
np <- import("numpy")
theme_set(theme_bw(base_size = 16))

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

### NSC
preds <- np$load("data/preds_targets/preds_NSC_Test.npy")
targets <- np$load("data/preds_targets/targets_NSC_Test.npy")

df <- cbind(preds, targets); colnames(df) <- c("preds", "targets")
df <- as.data.frame(df)
df$density <- get_density(df$targets, df$preds, h=c(1,1), n = 100)

ggplot(df, aes(x=targets, y=preds)) +
  geom_point(aes(x=targets, y=preds, color = density), size = .8) +
  sm_statCorr(show_text = F, color = "#BFACE2", size = 0.5) + 
  xlim(0,8) +
  ylim(0,8) + 
  theme_bw() +
  scale_color_gradientn(limits=c(0,1.5), breaks=seq(0, 1.5, by=0.3), colours=c( "grey80", "#a168d5" ,"#ff8100")) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  xlab("") +
  ylab("")

### ESC
preds <- np$load("data/preds_targets/preds_ESC_Valid.npy")
targets <- np$load("data/preds_targets/targets_ESC_Valid.npy")

df <- cbind(preds, targets); colnames(df) <- c("preds", "targets")
df <- as.data.frame(df)
df$density <- get_density(df$targets, df$preds, h=c(1,1), n = 500)

ggplot(df, aes(x=targets, y=preds)) +
  geom_point(aes(x=targets, y=preds, color = density), size = .8) +
  sm_statCorr(show_text = F, color = "#BFACE2", size = 0.5) + 
  xlim(0,6) +
  ylim(0,6) + 
  theme_bw() +
  scale_color_gradientn(limits=c(0,1.2), breaks=seq(0, 1.2, by=0.3), colours=c( "grey80", "#a168d5" ,"#ff8100")) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  xlab("") +
  ylab("")  

rm(list=ls())
gc()



##################
#####Fig. S5A#####
##################

### Fig. S5A: gnomAD rare variants enrichment

burden_dat <- read_excel("data/Supplementary Table 13 - gnomAD variants enrichment.xlsx", sheet = 1)
unique_groups <- unique(burden_dat$group)

subsets <- lapply(unique_groups, function(group_name) {
  subset(burden_dat, group == group_name)
})

results_list <- list()

for (i in seq_along(subsets)) {
  data <- subsets[[i]]
  
  # Check if the subset has only one row
  if (nrow(data) <= 1) {
    next
  }
  
  # Calculate odds ratio and 95% CI for each category compared to genome-wide
  odds_ratios <- vector("numeric", length = nrow(data) - 1)
  lower_ci <- vector("numeric", length = nrow(data) - 1)
  upper_ci <- vector("numeric", length = nrow(data) - 1)
  
  for (j in 2:nrow(data)) {
    # Create a 2x2 table
    variants <- matrix(c(data$rare_variants[j], data$rare_variants[1], 
                         data$size[j], data$size[1]), nrow = 2)
    
    # Calculate the odds ratio
    odds_ratio <- variants[1, 1] * variants[2, 2] / (variants[1, 2] * variants[2, 1])
    
    # Calculate standard error
    se_log_or <- sqrt(1 / variants[1, 1] + 1 / variants[2, 1] + 1 / variants[1, 2] + 1 / variants[2, 2])
    
    # Calculate confidence interval
    z_value <- qnorm(0.975)  # 95% confidence interval
    conf_interval <- exp(log(odds_ratio) + c(-1, 1) * z_value * se_log_or)
    
    # Store results
    odds_ratios[j - 1] <- odds_ratio
    lower_ci[j - 1] <- conf_interval[1]
    upper_ci[j - 1] <- conf_interval[2]
  }
  
  # Combine results
  result <- data.frame(
    index = 1:(nrow(data) - 1),
    category = data$region[-1],
    odds_ratio = odds_ratios,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    group = unique_groups[i]  # Add group column
  )
  
  # Store result in the results list
  results_list[[i]] <- result
}

combined_results <- do.call(rbind, results_list)

diagram <- c(c(paste0("chr",1:22), "chrX"),
             "codingExon", "random", "threeUTR", "fiveUTR", "DAE", "nDAE", "ENCODE_cCREs", "VISTA", 
             "common_high", "ESC_high", "NSC_high",
             "ESC_category_1", "ESC_category_2", "ESC_category_3", 
             "ESC_category_4", "ESC_category_5", "ESC_top10",
             "NSC_category_1", "NSC_category_2", "NSC_category_3", 
             "NSC_category_4", "NSC_category_5", "NSC_top10") # all in genomAD

plot_dat <- subset(combined_results, group == "gnomAD" & category %in% diagram) 
gnomAD_OE <- plot_dat[plot_dat$group == "gnomAD", "odds_ratio"]

plot_dat$category <- factor(plot_dat$category, levels = diagram)
ggplot(plot_dat, aes(y = category, x = odds_ratio)) +
  geom_point(size = 2, color = "#0E4c92") +  
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0, color = "#0E4c92") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
  # scale_y_continuous(name = "", breaks=1:nrow(plot_dat), labels = plot_dat$category, trans = "reverse") +
  xlim(c(0.70,1.25)) + 
  xlab("Odds Ratio (95% CI)") + 
  ylab("") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.text.x.bottom = element_text(size = 8, colour = "black"),
        axis.title.x = element_text(size = 8, colour = "black")) +
  geom_text(aes(label = sprintf("%.6f", odds_ratio)), hjust = 0.6, vjust = -0.7, size = 2.5)

rm(list=ls())
gc()
