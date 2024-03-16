library(dplyr)
library(ggplot2)
# library(ggpubr)
library(ggrepel)
library(readxl)
library(tidyr)
library(data.table)
library(readxl)
options(scipen = 999)
setwd("/Users/Ruizhi/Work/EMC/Projects/Deep_learning/motifs/")

##########################
###Motifs visualization###
##########################


motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")

# NSC and ESC need to be loaded
tomtom <- read.table("../data/patterns/NSC_ESC_common/tomtom_results_ESC_high_NSCmodel.txt",sep = "\t") %>% .[,1:6]; colnames(tomtom) <- c("pattern", "num", "Motif", "pval", "eval", "qval")
tomtom <- merge(motif2gene, tomtom, by = "Motif")
tomtom <- tomtom %>%
  mutate(nor_eval = -log10(eval)) %>% 
  mutate(nor_pval = -log10(pval)) 
tomtom$pattern <- gsub("\\..*", "", tomtom$pattern)

ggplot(tomtom, aes(x = num, y = nor_pval)) + 
  geom_point(aes(color = pattern, size = num)) +
  scale_size(range = c(3, 3))  + ggrepel::geom_text_repel(aes(label=GeneName), size=4, max.overlaps = 25) +
  xlab("The number of seqlets") + ylab("-log10(pvalue)")  + theme_bw() +
  theme(legend.position = "none") +
  xlim(0, 920) +
  ylim(2, 12.5)  +
  scale_color_manual(values=c("#234471", "#cb1517")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        panel.grid=element_blank())

#####################
###PRROC for AUPRC###
#####################
library(PRROC)
preds_targets<- read.csv("../data/Top_Bottom_ESC.csv")
# calculate precision and recall values
pr <- pr.curve(scores.class0 = preds_targets$preds, weights.class0 = preds_targets$label, curve = TRUE)
plot(pr)

roc<-roc.curve(scores.class0 = preds_targets$preds, weights.class0 = preds_targets$label, curve = TRUE)
plot(roc)

#####################
###PRROC for AUPRC###
#####################
library(smplot2)
library(ggrastr)
library(scales)
edi_high <- read.csv("../data/preds_targets_ediESC.csv")
edi_high$re_targets <- rescale(edi_high$targets, to = c(min(edi_high$preds), max(edi_high$preds)))

ggplot(data = edi_high, mapping = aes(x = re_targets, y = preds)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
  sm_statCorr() +
  ggrastr::rasterise(geom_point(size=0.1, alpha=0.1))


################################
###preds, targets correlation###
################################
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

# NSC
preds <- np$load("../data/preds_targets/preds_NSC_Test.npy")
targets <- np$load("../data/preds_targets/targets_NSC_Test.npy")

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
  
# ESC
preds <- np$load("../data/preds_targets/preds_ESC_Valid.npy")
targets <- np$load("../data/preds_targets/targets_ESC_Valid.npy")

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
  
##############################
###contribution scores plot###
##############################
library(ggseqlogo)

# NSC
hypoth_scores <- np$load("../data/contri_scores/shap_explanations_NSC.npy")
inps <- np$load("../data/contri_scores/inp_NSC.npy")

contri_scores <- np$multiply(hypoth_scores, inps)

# one example
examp <- contri_scores[13776 + 1,1:4, 1:1000]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 1000, by = 100), expand = c(0, 0), limits = c(0, NA)) + 
  annotate('rect', xmin = 471, xmax = 482, ymin = 0, ymax = 0.03, alpha = .1, fill='#ff3399', size = 0.3) +
  annotate('rect', xmin = 579, xmax = 596, ymin = 0, ymax = 0.03, alpha = .1, fill='#ff3399', size = 0.3) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none")

# zoom in, 7.5 x 1.2
examp <- contri_scores[13776 + 1,1:4, 470:600]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 130, by = 10), expand = c(0, 0), limits = c(0, NA)) + 
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none")



# ESC
hypoth_scores <- np$load("../data/contri_scores/shap_explanations_ESC.npy")
inps <- np$load("../data/contri_scores/inp_ESC.npy")

contri_scores <- np$multiply(hypoth_scores, inps)

# one example, padded
examp <- contri_scores[6511 + 1,1:4, 1:1000]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) + 
  annotate('rect', xmin = 503, xmax = 518, ymin = -0.03, ymax = 0.01, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  annotate('rect', xmin = 617, xmax = 630, ymin = -0.03, ymax = 0.01, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))

# one example, unpadded
examp <- contri_scores[6511 + 1,1:4, 199:801]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) + 
  annotate('rect', xmin = 305, xmax = 320, ymin = -0.03, ymax = 0.01, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  annotate('rect', xmin = 419, xmax = 432, ymin = -0.03, ymax = 0.01, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))

examp <- contri_scores[6511 + 1,1:4, 502:632]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 130, by = 10), expand = c(0, 0), limits = c(0, NA)) + 
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none")


###############################
###############################
###contribution scores plot#### I have to repeat it again!!!!!! the NSC_ESC_common shit!
###############################
###############################
library(ggseqlogo)

# NSC
hypoth_scores <- np$load("../data/contri_scores/NSC_ESC_common/shap_explanations_NSC_high.npy")
inps <- np$load("../data/contri_scores/NSC_ESC_common/inp_NSC_high.npy")

contri_scores <- np$multiply(hypoth_scores, inps)

# one example
examp <- contri_scores[7375 + 1,1:4, 1:1000]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 1000, by = 100), expand = c(0, 0), limits = c(0, NA)) + 
  annotate('rect', xmin = 465, xmax = 486, ymin = 0, ymax = 0.03, alpha = .1, fill='#ff3399', size = 0.3) +
  annotate('rect', xmin = 579, xmax = 597, ymin = 0, ymax = 0.03, alpha = .1, fill='#ff3399', size = 0.3) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none")

# zoom in, 7.5 x 1.2
examp <- contri_scores[7375 + 1,1:4, 465:600]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 130, by = 10), expand = c(0, 0), limits = c(0, NA)) + 
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none")



# ESC
hypoth_scores <- np$load("../data/contri_scores/NSC_ESC_common/shap_explanations_ESC_high.npy")
inps <- np$load("../data/contri_scores/NSC_ESC_common/inp_ESC_high.npy")

contri_scores <- np$multiply(hypoth_scores, inps)

# one example, padded
examp <- contri_scores[3451 + 1,1:4, 216:784]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) + 
  annotate('rect', xmin = 63, xmax = 93, ymin = -0.01, ymax = 0.02, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  annotate('rect', xmin = 238, xmax = 268, ymin = -0.01, ymax = 0.02, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))

# one example, unpadded
examp <- contri_scores[3451 + 1,1:4, 199:801]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 1000, by = 100)) + 
  annotate('rect', xmin = 305, xmax = 320, ymin = -0.03, ymax = 0.01, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  annotate('rect', xmin = 419, xmax = 432, ymin = -0.03, ymax = 0.01, alpha = .1, col='black', fill='#ff3399', size = 0.3) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))

examp <- contri_scores[6511 + 1,1:4, 502:632]; rownames(examp) <- c('A', 'C', 'G', 'T')
ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
  scale_x_continuous(breaks = seq(0, 130, by = 10), expand = c(0, 0), limits = c(0, NA)) + 
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none")


##########################
###comprehensive motifs###
##########################


######################
###examples of ncER###
######################

#ESC
examp_ncER <- read.table("../motifs/example_ncER.txt", sep = "\t", header = F) %>%
  subset(V1 %in% "chr2") 
colnames(examp_ncER) <- c("chr","start", "end", "ncER")
examp_ncER$group <- rep(1:100, each = 15)[1:nrow(examp_ncER)]

# ggplot(examp_ncER, aes(end, ncER)) + 
#   geom_line()+ 
#   scale_x_continuous(breaks = seq(122247643, 122248245, by = 60), expand = c(0, 0)) +
#   ylim(0,100) +
#   geom_vline(xintercept=305+122247643, linetype="dashed", 
#              color = "red", size=.2) +
#   geom_vline(xintercept=320+122247643, linetype="dashed", 
#              color = "red", size=.2) +
#   geom_vline(xintercept=419+122247643, linetype="dashed", 
#              color = "red", size=.2) +
#   geom_vline(xintercept=432+122247643, linetype="dashed", 
#              color = "red", size=.2) +
#   theme_classic() 

ggplot(examp_ncER, aes(x=group, y=ncER, group = group)) + 
  geom_boxplot(outlier.size=0.2, fill="#3399ff", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  ylim(0,100)

# NSC
examp_ncER <- read.table("../motifs/example_ncER.txt", sep = "\t", header = F) %>%
  subset(V1 %in% "chr8") 

colnames(examp_ncER) <- c("chr","start", "end", "ncER")
examp_ncER$group <- rep(1:100, each = 25)[1:nrow(examp_ncER)]

# ggplot(examp_ncER, aes(start, ncER)) + 
#   geom_line()+ 
#   scale_x_continuous(breaks = seq(124491517, 124492516, by = 100), expand = c(0, 0)) +
#   ylim(0,100) +
#   theme_classic() +
#   geom_vline(xintercept=503+124491517, linetype="dashed", 
#              color = "red", size=.2) +
#   geom_vline(xintercept=518+124491517, linetype="dashed", 
#              color = "red", size=.2) +
#   geom_vline(xintercept=617+124491517, linetype="dashed", 
#            color = "red", size=.2) +
#   geom_vline(xintercept=630+124491517, linetype="dashed", 
#              color = "red", size=.2)

ggplot(examp_ncER, aes(x=group, y=ncER, group = group)) + 
  geom_boxplot(outlier.size=0.2, fill="#ff3399", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=12),
        legend.position="none") +
  ylim(0,100)
  
############################
###Motifs gene expression###
############################
library(tidyverse)
library(gplots)
library(RColorBrewer)
# load wt gene expression of NSC and ESC 
wt_expression <- read.csv("../data/wt_Expression.RPKM.csv") %>% 
  select(c("GeneID", "GeneName", "ESC.YY1_mean", "NSC.YY1_mean"))
colnames(wt_expression) <- c("GeneID","GeneName", "ESC", "NSC")

# load tomtom result of NSC and ESC 
motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")

tomtom <- read.table("../data/patterns/NSC_ESC_common/tomtom_results_ESC_high.txt",sep = "\t") %>% .[,1:6]; colnames(tomtom) <- c("pattern", "num", "Motif", "pval", "eval", "qval")
tomtom <- merge(motif2gene, tomtom, by = "Motif")
tomtom <- tomtom %>%
  mutate(nor_eval = -log10(eval)) %>% 
  mutate(nor_pval = -log10(pval)) 
tomtom$pattern <- gsub("\\..*", "", tomtom$pattern)

tomtom$GeneName <- toupper(tomtom$GeneName)
tomtom <- separate_rows(tomtom, GeneName, sep = "::")
tomtom <- separate_rows(tomtom, GeneName, sep = "-")

# motif expression plot
tomtom_expr <- tomtom %>% 
  left_join(wt_expression, by = "GeneName") %>%
  select(c("GeneName", "pattern", "num", "ESC", "NSC")) %>%
  arrange(desc(pattern), desc(num)) %>%
  select(c("GeneName", "ESC", "NSC")) %>%
  unique() %>% as.data.frame()
rownames(tomtom_expr) <-tomtom_expr$GeneName
tomtom_expr <- tomtom_expr[-1]
tomtom_expr <- as.matrix(log2(tomtom_expr + 1))

# plot
# heatmap.2(tomtom_expr,  na.color = "white", trace ="none", 
#           Rowv=FALSE, Colv=FALSE, dendrogram = "none", cexRow=0.5,cexCol=1, margins=c(6,6),
#           col = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(20))

# colfunc <- colorRampPalette(c("#5690cc", "#74a4d5", "#92b8de", "#b0cbe7", "#cfdff0" ,"#edf3f9",
#                               "#fef2f8", "#fbcee4", "#f8aad1", "#f486bd", "#f162aa", "#ee3e96"))
colfunc <- colorRampPalette(c("#234471", "#0962A6", "#0680C3", "#46A8DF", "#87C6EE" ,"#C4E5F7",
                              "#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#FF0000", "#CC0000"))
Breaks <- seq(0,8.2, length = 21)
heatmap.2(tomtom_expr,  na.color = "white", trace ="none", 
          Rowv=FALSE, Colv=FALSE, dendrogram = "none", cexRow=0.5,cexCol=1, margins=c(6,6),
          col = colfunc(20), breaks = Breaks)



############################
###patterns linked genes####
############################
library("readxl")

######
###### NSC
######
# load tomtom result
motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")
tomtom <- read.table("tomtom_results_High_NSC.txt",sep = "\t") %>% .[,1:6]; colnames(tomtom) <- c("pattern", "num", "Motif", "pval", "eval", "qval")
tomtom <- merge(motif2gene, tomtom, by = "Motif")
tomtom <- separate(tomtom, col = pattern, into = c("pattern_type", "pattern_name"), sep = "\\.")
tomtom<-arrange(tomtom, desc(pattern_type), pattern_name)
# load patterns' coordinates
patterns <- read.table("../data/patterns/patterns_coordinates_NSC.txt", sep = "\t", header = T)
patterns <- left_join(patterns, tomtom[,1:4], by = c("pattern_type", "pattern_name"))
colnames(patterns) <- c("pattern_type","pattern_name","enhancers","Motif","Motif_GeneName")    

# load wt gene expression talbe
wt_expression <- read.csv("../data/wt_Expression.RPKM.csv") %>% 
  select(c("GeneID", "GeneName", "ESC.YY1_mean", "NSC.YY1_mean"))
colnames(wt_expression) <- c("GeneID","GeneName", "ESC", "NSC")

# load enhancers' linked genes
enhancer_genes <- read_excel("../data/patterns/Scaffold_Categories_TargetGenes.xlsx", sheet = 1)
patterns_genes <- patterns %>% left_join(enhancer_genes[,c(1,6:8)], by = "enhancers") %>%
  left_join(wt_expression, by = c("GeneID","GeneName")) %>%
  subset(GeneType %in% "protein-coding") 

patterns_plot <- patterns_genes %>%
  select(c("pattern_type", "pattern_name", "enhancers", "GeneName", "ESC", "NSC" )) %>%
  unique() %>%
  gather(key = "Cell_Type", value = "RPKM", -pattern_type, -enhancers, -GeneName, -pattern_name) %>%
  mutate(trans_RPKM=log2(RPKM + 1)) 

patterns_plot$pattern_name <- factor(patterns_plot$pattern_name, levels = paste0("pattern_",0:22))


ggplot(patterns_plot, aes(x=pattern_name, y=trans_RPKM, fill=Cell_Type)) + 
  geom_boxplot(outlier.size=0.2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~pattern_type) +
  ylim(0,10) +
  xlab("") +
  ylab("log2(RPKM + 1)") +
  scale_fill_manual(values = c("#5690CC", "#EE3E96"))



######
###### ESC
######
# load tomtom result
motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")
tomtom <- read.table("tomtom_results_High_ESC.txt",sep = "\t") %>% .[,1:6]; colnames(tomtom) <- c("pattern", "num", "Motif", "pval", "eval", "qval")
tomtom <- merge(motif2gene, tomtom, by = "Motif")
tomtom <- separate(tomtom, col = pattern, into = c("pattern_type", "pattern_name"), sep = "\\.")

# load patterns' coordinates
patterns <- read.table("../data/patterns/patterns_coordinates_ESC.txt", sep = "\t", header = T)
patterns <- left_join(patterns, tomtom[,1:4], by = c("pattern_type", "pattern_name"))
colnames(patterns) <- c("pattern_type","pattern_name","enhancers","Motif","Motif_GeneName")    

# load wt gene expression talbe
wt_expression <- read.csv("../data/wt_Expression.RPKM.csv") %>% 
  select(c("GeneID", "GeneName", "ESC.YY1_mean", "NSC.YY1_mean"))
colnames(wt_expression) <- c("GeneID","GeneName", "ESC", "NSC")

# load enhancers' linked genes
enhancer_genes <- read_excel("../data/patterns/Scaffold_Categories_TargetGenes.xlsx", sheet = 2)
patterns_genes <- patterns %>% left_join(enhancer_genes[,c(1,6:8)], by = "enhancers") %>%
  left_join(wt_expression, by = c("GeneID","GeneName")) %>%
  subset(GeneType %in% "protein-coding") 

patterns_plot <- patterns_genes %>%
  select(c("pattern_type", "pattern_name", "enhancers", "GeneName", "ESC", "NSC" )) %>%
  unique() %>%
  gather(key = "Cell_Type", value = "RPKM", -pattern_type, -enhancers, -GeneName, -pattern_name) %>%
  mutate(trans_RPKM=log2(RPKM + 1)) 

patterns_plot$pattern_name <- factor(patterns_plot$pattern_name, levels = paste0("pattern_",0:28))

ggplot(patterns_plot, aes(x=pattern_name, y=trans_RPKM, fill=Cell_Type)) + 
  geom_boxplot(outlier.size=0.2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(~pattern_type)+
  ylim(0,10) +
  xlab("") +
  ylab("log2(RPKM + 1)") +
  scale_fill_manual(values = c("#5690CC", "#EE3E96"))



#########################################
###prepare regions for epigenome data####
#########################################

# only separate NSC/ESC
patterns <- read.table("../data/patterns/patterns_coordinates_NSC.txt", sep = "\t", header = T)
regions <- patterns %>% 
  separate(col = location, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  select(c("Chr", "Start", "End")) %>%
  unique() %>% 
  arrange(Chr, as.numeric(Start))

write.table(regions, "../data/patterns/Coordinates_ESC_High.bed", row.names = F, quote = F,col.names = F, sep = "\t")  

# separate NSC neg + pos, ESC neg
patterns <- read.table("../data/patterns/patterns_coordinates_NSC.txt", sep = "\t", header = T)
regions <- patterns %>% 
  subset(pattern_type == "neg_patterns") %>% # neg_patterns, pos_patterns
  separate(col = location, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  select(c("Chr", "Start", "End")) %>%
  unique() %>% 
  arrange(Chr, as.numeric(Start))

write.table(regions, "../data/patterns/Coordinates_NSC_High.Neg.bed", row.names = F, quote = F,col.names = F, sep = "\t") 

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


############################################
###add HPO information to enhancer tables###
############################################
# add HPO information to enhancer tables

# enhancer_category excel
excel_file <- "/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets.xlsx"
sheets <- excel_sheets(excel_file)
enhancer_category <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(enhancer_category) <- sheets

# Generate the complete OMIN_HPO_GeneName table
omim <- read.delim("../data/HPO_OMIM/OMIM_GeneMap_21.04.2022.txt", sep="\t") %>%
  select(c("Approved.Gene.Symbol","Phenotypes")) %>%
  subset(.$Phenotypes != "")  %>%
  subset(.$Approved.Gene.Symbol != "")

hpo_to_genes <- read.table("../data/HPO_OMIM/phenotype_to_genes.txt", sep="\t", header = T) %>% 
  dplyr::select(c("gene_symbol","hpo_id")) %>%
  unique()

omim_hpo_gene <- omim %>% left_join(hpo_to_genes, by=c('Approved.Gene.Symbol'='gene_symbol'), relationship = "many-to-many") %>%
  group_by(Approved.Gene.Symbol, Phenotypes) %>%
  summarize(hpo_id = paste(hpo_id, collapse = ";")) #summarize(hpo_id = toString(hpo_id))
colnames(omim_hpo_gene) <- c("GeneName", "OMIM_Phenotypes", "HPO" )

# left_join the enhancer_category with hpo information
enhancer_category_hpo <- enhancer_category %>%
  lapply( function(x) dplyr::left_join(x, omim_hpo_gene, by="GeneName")) 

# save as excel enhancer_category_hpo
library(openxlsx)
excel_file <- "/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx"
wb <- createWorkbook()

for (i in seq_along(enhancer_category_hpo)) {
  sheet_name <- names(enhancer_category_hpo)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, enhancer_category_hpo[[i]])
}

saveWorkbook(wb, excel_file, overwrite = T)

######################################################
###prepare motif location talbes with targeted gene###
######################################################

NSC_whole_table <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>% 
  .[,c(1,8:10)]; colnames(NSC_whole_table) <- c("location", "NSC_TargetGene", "NSC_OMIM_Phenotypes", "NSC_HPO")

ESC_whole_table <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>%
  .[,c(1,8:10)]; colnames(ESC_whole_table) <- c("location", "ESC_TargetGene", "ESC_OMIM_Phenotypes", "ESC_HPO")

# enhancer_motif excel
excel_file <- "../data/patterns/NSC_ESC_common/enhancer_motif.xlsx"
sheets <- excel_sheets(excel_file)
enhancer_motif <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(enhancer_motif) <- sheets

# tomtom_result excel
excel_file <- "../data/patterns/NSC_ESC_common/tomtom_result.xlsx"
sheets <- excel_sheets(excel_file)
tomtom_result <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(tomtom_result) <- sheets

# motif -> GeneName table
motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")

# convert MotifName to GeneName
tomtom_GeneName <- tomtom_result %>% lapply( function(x) tidyr::separate(x, col = pattern, into = c("pattern_type", "pattern_name"), sep = "\\.")) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match0" = "Motif"))) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match1" = "Motif"))) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match2" = "Motif"))) %>%
  lapply( function(x) dplyr::select(x, c("pattern_type", "pattern_name", "num_seqlets", "GeneName.x", "GeneName.y", "GeneName"))) %>%
  lapply( function(x) setNames(x, c("pattern_type", "pattern_name", "num_seqlets", "motif_GeneName0", "motif_GeneName1", "motif_GeneName2")))
  
# left_join enhancer_motif with Motif's GeneName  
enhancer_motifGeneName <- Map(function(x, y) dplyr::left_join(x, y, by=c("pattern_type", "pattern_name")), enhancer_motif, tomtom_GeneName)
  
enhancer_motifGeneName <- enhancer_motifGeneName %>% lapply( function(x) dplyr::left_join(x, NSC_whole_table, "location", relationship = "many-to-many" )) %>%
  lapply( function(x) dplyr::left_join(x, ESC_whole_table, "location", relationship = "many-to-many" )) %>%
  lapply( function(x) x[,c(1,9:14,2:8)])

# save as excel
library(openxlsx)
excel_file <- "../data/patterns/NSC_ESC_common/enhancer_motif_GeneName.xlsx"
wb <- createWorkbook()

for (i in seq_along(enhancer_motifGeneName)) {
  sheet_name <- names(enhancer_motifGeneName)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, enhancer_motifGeneName[[i]])
}

saveWorkbook(wb, excel_file, overwrite = T)


###################################################
###Prepare selected_enhancers for lab validation###
###################################################
wt_expression <- read.csv("../data/wt_Expression.RPKM.csv") %>% 
  select(c("GeneName", "NSC.YY1_mean"))
colnames(wt_expression) <- c("GeneName", "NSC_RPKM")


selected_HPO <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Supplementary Tables/HPO_terms_selected_final.xlsx", col_types = c(rep("text", 3)))

NSC_high <- read_excel("../data/patterns/NSC_ESC_common/enhancer_motif_GeneName.xlsx", sheet="NSC_high", col_types = c(rep("text", 14))) %>%
  dplyr::filter(!is.na(pattern_name)) %>%
  subset(NSC_TargetGene %in% selected_HPO$GeneName)  %>%
  subset(pattern_type == "pos_patterns") %>%
  subset(motif_GeneName0 %in% c("TP53", "ZFP42")) %>%
  .[,-c(5,6,7)] %>% unique() %>%
  left_join(wt_expression, by = c("NSC_TargetGene" = "GeneName")) %>%
  .[,c(1:2,12,3:11)]
write.csv(NSC_high, "../data/patterns/NSC_ESC_common/NSC_high_selectedCoordinates.csv", row.names = F)
  
exclude_enhancer <- read.table("../data/patterns/NSC_ESC_common/exclude_enhancer.txt", header = T)
NSC_high <- NSC_high %>%
  subset(! location %in% exclude_enhancer$location)
write.csv(NSC_high, "../data/patterns/NSC_ESC_common/NSC_high_selectedCoordinates.csv", row.names = F)



###########################################
###split enhancer activity to percentile###
###########################################
# activity <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets.xlsx") %>%
#   .[,-c(6:8)] %>% unique() %>%
#   mutate(percentile_NSC = (rank(`log2(avg.NSC.RPP+1)`) / length(`log2(avg.NSC.RPP+1)`))*100) %>%
#   mutate(percentile_ESC = (rank(`log2(avg.ESC.RPP+1)`) / length(`log2(avg.ESC.RPP+1)`))*100) %>%
#   .[,c(1,6,7)] %>%
#   separate(col = enhancers, into = c("Chr", "Start_tmp"), sep = ":") %>%
#   separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
#   arrange(Chr, as.numeric(Start))

activity <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets.xlsx") %>%
  .[,-c(6:8)] %>% unique() %>% 
  arrange(as.numeric(`log2(avg.NSC.RPP+1)`)) %>%
  mutate(percentile_NSC = 1:nrow(.) )  %>%
  arrange(as.numeric(`log2(avg.ESC.RPP+1)`)) %>%
  mutate(percentile_ESC = 1:nrow(.) ) %>%
  .[,c(1,6,7)] %>%
  separate(col = enhancers, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  arrange(Chr, as.numeric(Start))

write.table(activity, "../data/Percentile_approach/Enhancer_percent.txt", sep = "\t", quote = F, row.names = F)


#####################################
###plot percentile motifs patterns###
#####################################
options(scipen = 999)
# tomtom_result excel
excel_file <- "../data/Percentile_approach/tomtom_result_activity.xlsx"
# excel_file <- "../data/patterns/NSC_ESC_common/tomtom_result.xlsx"
sheets <- excel_sheets(excel_file)
col_types <- c("text", "numeric", rep(c("text", "numeric", "numeric", "numeric"), 3))
tomtom_result <- lapply(sheets, function(x) read_excel(excel_file, sheet = x, col_types = col_types))
names(tomtom_result) <- sheets

# motif -> GeneName table
motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")

# convert MotifName to GeneName
tomtom_GeneName <- tomtom_result %>% 
  lapply( function(x) x[,1:6]) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match0" = "Motif"))) %>%
  lapply( function(x) mutate(x, nor_eval = -log10(eval0) )) %>%
  lapply( function(x) mutate(x, nor_pval = -log10(pval0) )) %>%
  lapply( function(x) mutate(x, pattern = gsub("\\..*", "", x$pattern) )) 

plots <- lapply(names(tomtom_GeneName), function(name) {
  ggplot_obj <- tomtom_GeneName[[name]]
  
  ggplot_obj$font_face <- ifelse(ggplot_obj$qval0 < 0.05, "bold", "plain")
  
  ggplot(ggplot_obj, aes(x = num_seqlets, y = nor_pval)) + 
    geom_point(aes(color = pattern, size = num_seqlets)) +
    scale_size(range = c(3, 3)) +
    ggrepel::geom_text_repel(aes(label = GeneName, fontface = font_face), size = 4, max.overlaps = 35) +
    xlab("The number of seqlets") + ylab("-log10(pvalue)") +
    theme_bw() +
    theme(legend.position = "none") +
    xlim(0, 1250) +
    ylim(2, 13) +
    scale_color_manual(values = c("#234471", "#cb1517")) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          panel.grid = element_blank()) +
    ggtitle(name)  # Add the list's name as the title
})


pdf('~/Downloads/percentile_motif_activity.pdf', onefile = TRUE)
# pdf('~/Downloads/activity_motif.pdf', onefile = TRUE)
for (i in seq(length(plots))) {
  print(plots[[i]])
}
dev.off()

#################################
###plot motif polygenetic tree###
#################################
library(readxl)
library(dplyr)
library(reshape2)
library(MotifDb)
library(motifStack)

# generate the selected motifs from subgroups
tomtom <- read_excel("../data/patterns/NSC_ESC_common/tomtom_result.xlsx") %>% select("match0", "match1", "match2") %>% .[1:10,]
tomtom_l <- unlist(tomtom, use.names = FALSE)
tomtom_l <- tomtom_l[!is.na(tomtom_l)]
tomtom_l <- lapply(tomtom_l, function(x) paste(x, ".jaspar", sep = ""))

tomtom_l <- c("MA")

# read pcm from the selected motifs in JASPAR
path <- "../data/JASPAR2022_CORE_non-redundant_pfms_jaspar/"
files <- list.files(path, pattern = "\\.jaspar$", full.names = TRUE)
files <- files[basename(files) %in% tomtom_l]
motifs <- importMatrix(files, format = "pcm", to = "pfm")

motifPiles(motifs)

motifStack(motifs, layout="phylog", f.phylog=.15, f.logo=0.35)



############################
###add GeneName to Motifs###
############################
library(dplyr)
library(readxl)
# tomtom_result excel
excel_file <- "Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary_Table_13 - Comprehensive motifs identified from TF-modisco-lite.xlsx"
sheets <- excel_sheets(excel_file)
tomtom_result <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(tomtom_result) <- sheets

# motif -> GeneName table
motif2gene <- read.table("Work/EMC/Projects/Deep_learning/motifs/motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")

# convert MotifName to GeneName
tomtom_GeneName <- tomtom_result %>% 
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match0" = "Motif"))) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match1" = "Motif"))) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match2" = "Motif"))) %>%
  lapply( function(x) x[,c(1:3,15,4:7,16,8:11,17)]) %>%
  lapply( function(x) setNames(x, c("pattern", "num_seqlets", "match0", "match0_GeneName", "pval0", "eval0", "qval0", 
                                    "match1", "match1_GeneName", "pval1", "eval1",  "qval1", "match2", "match2_GeneName")))
# save as excel
library(openxlsx)
excel_file <- "Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 13 - Comprehensive motifs identified from TF-modisco-lite.xlsx"
wb <- createWorkbook()

for (i in seq_along(tomtom_GeneName)) {
  sheet_name <- names(tomtom_GeneName)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tomtom_GeneName[[i]])
}

saveWorkbook(wb, excel_file, overwrite = T)



#########################################################
###correlatioin of ChIP-STARR, prediction and the FACS###
#########################################################
library(reshape2)
library(ggpubr)
# read the FACS of validation enhancers
facs <- read_excel("../data/Experiment/FACS-NSC-GFP.xlsx") %>% 
  t() %>% as.data.frame() %>% 
  mutate(FACS = rowMeans(.) * 0.01) %>%
  mutate(NSC_TargetGene = row.names(.)) %>%
  dplyr::select("NSC_TargetGene", "FACS")

# ChIP-STARR-seq and Prediction
selected_enhancer <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Version1/Supplementary Tables/Supplementary Table 15 - Selected validated enhancers.xlsx", sheet = 2)
whole_nsc <- read.csv("/Users/Ruizhi/Downloads/whole_NSC.csv")

symbol <- c("IRF2BPL", "DPYSL5", "KNL1", "ADAR", "TRIO", "OAT", "DNMT3A", "CIC","TUBB2B", "ACTB", "ZMIZ1", "KDM4B", "NAA20", "NAT8L", "ZBTB11", "LAMB2", "TKT", "ASH1L", "PAFAH1B1", "ATP6V1A")
selected_enhancer <- selected_enhancer %>% left_join(whole_nsc, by = "location") %>%
  subset(NSC_TargetGene %in% symbol) %>% 
  left_join(facs, by = "NSC_TargetGene") %>%
  .[,c(1,2,13:15)]

selected_enhancer$NSC_TargetGene <- factor(selected_enhancer$NSC_TargetGene, levels = symbol)

# plot figures
library(tidyr)
activity <- gather(selected_enhancer, key = "group", value = "activity", -c("location", "NSC_TargetGene"))
activity$group <- factor(activity$group, levels = c("targets", "preds", "FACS"))
bar_activity <- ggplot(activity, aes(x = NSC_TargetGene, y = activity, fill=group)) + 
  geom_bar(stat="identity", color="black", position=position_dodge())+
  ylab("") +
  xlab("") + 
  scale_fill_discrete(labels=c('ChIP-STARR', 'Prediction', 'FACS')) + 
  theme(axis.title = element_text(size = 8),
        axis.text=element_text(size=8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0))+
  expand_limits(y = c(0,4))

corr_TargetFacs <- ggscatter(selected_enhancer[,3:5], 
                  x = "targets", y = "FACS",
                  add = "reg.line", 
                  conf.int = TRUE,
                  palette = "jco",cor.method = "pearson") +
  stat_cor() + 
  xlab("ChIP-STARR-seq") +
  ylab("FACS") +
  theme(axis.title = element_text(size = 10),
        axis.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0))

corr_TargetPred <- ggscatter(selected_enhancer[,3:5], 
                             x = "preds", y = "FACS",
                             add = "reg.line", 
                             conf.int = TRUE,
                             palette = "jco",cor.method = "pearson") +
  stat_cor() + 
  xlab("Prediction") +
  ylab("FACS") +
  theme(axis.title = element_text(size = 10),
        axis.text=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0))

ggarrange(bar_activity,                                                 # First row with scatter plot
          ggarrange(corr_TargetFacs, corr_TargetPred, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 


#######################################################
###contribution scores plots of validation enhancers###
#######################################################

#### contribution scores
library(ggseqlogo)
# read the enhancer location which match the order of contribution score
enhancer_activity <- read.table("../data/contri_scores/Enhancer_activity.txt",sep = "\t", header = T) %>%
  unite("tmp",Chr:Start, sep = ":") %>%
  unite("location", tmp:End, sep = "-")
# NSC
hypoth_scores <- np$load("../data/contri_scores/shap_explanations_NSC.npy")
inps <- np$load("../data/contri_scores/inp_NSC.npy")

contri_scores <- np$multiply(hypoth_scores, inps)

# Find the row name (ID) where the value matches
selected_enhancer <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Version1/Supplementary Tables/Supplementary Table 15 - Selected validated enhancers.xlsx", sheet = 3)

enhancer_ID <- selected_enhancer$location
# enhancer_order <- as.numeric(row.names(enhancer_activity)[which(enhancer_activity$location == enhancer_ID)])
# enhancer_order <- rownames(enhancer_activity[enhancer_activity$location %in% enhancer_ID, ])

# one example
library(gridExtra)

plots <- list()


for (i in 1:length(enhancer_ID) ) {
  
  enhancer_order <- as.numeric(row.names(enhancer_activity)[which(enhancer_activity$location == enhancer_ID[i])])
  
  examp <- contri_scores[enhancer_order,1:4, 1:1000]; rownames(examp) <- c('A', 'C', 'G', 'T')
  
  plots[[i]] <- ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
    scale_x_continuous(breaks = seq(0, 1000, by = 100), expand = c(0, 0), limits = c(0, NA)) + 
    ggtitle(enhancer_ID[i]) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.5),
          axis.text=element_text(size=8,colour = "black"),
          legend.position="none",
          plot.title = element_text(color="black", size=8))
  
}
pdf("../figures/sequence_logos.pdf", width = 7.63, height = 22)

arranged_plots <- do.call(grid.arrange, c(plots[1:20], ncol = 1))

dev.off()



#### ncER scores

validation_ncER <- read.table("../data/contri_scores/Validate_enhancers.hg19.ncER.txt", sep = "\t", header = F)
colnames(validation_ncER) <- c("chr","start", "end", "ncER")

selected_enhancer_hg19 <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Version1/Supplementary Tables/Supplementary Table 15 - Selected validated enhancers.xlsx", sheet = 4)
enhancer_ID <- selected_enhancer_hg19$location
motif_ID <- selected_enhancer_hg19$motif_location

plots_ncER <- list()


for (i in 1:length(enhancer_ID) ) {
  # options(repr.plot.width =9, repr.plot.height =9)
  # return enhancer location
  split_enhancer <- strsplit(enhancer_ID[i], "[:\\-]")
  chr_enhanhcer <- split_enhancer[[1]][1]
  start_enhanhcer <- as.numeric(split_enhancer[[1]][2])
  end_enhanhcer <- as.numeric(split_enhancer[[1]][3])
  
  # return motif location
  split_motif <- strsplit(motif_ID[i], "[:\\-]")
  chr_motif <- split_motif[[1]][1]
  start_motif <- as.numeric(split_motif[[1]][2])
  end_motif <- as.numeric(split_motif[[1]][3])  
  
  ncER <- subset(validation_ncER, validation_ncER$chr == chr_enhanhcer & validation_ncER$start >= start_enhanhcer & validation_ncER$end <= end_enhanhcer)
  ncER_in_motif <- subset(validation_ncER, validation_ncER$start >= start_motif & validation_ncER$end <= end_motif)
  
  # ncER$group <- rep(1:100, each = 25)[1:nrow(ncER)]
  
  plots_ncER[[i]] <- ggplot(ncER, aes(x=start, y=ncER)) + 
    geom_line(size = 0.3)+
    geom_point(data = ncER_in_motif, aes(x = start, y = ncER), color = "red", size = 0.2) +  # Highlighted points in red
    scale_x_continuous(expand = c(0, 0)) +
    xlab("") +
    ggtitle(enhancer_ID[i]) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.5),
          axis.text=element_text(size=8,colour = "black"),
          legend.position="none",
          plot.title = element_text(color="black", size=8)) +
    ylim(0,100)
  
}


pdf("../figures/ncER_hg19.pdf", width = 7.63, height = 22)

arranged_plots_ncER <- do.call(grid.arrange, c(plots_ncER[1:20], ncol = 1))

dev.off()


pdf("../figures/combined_plots.pdf", width = 7.63, height = 44)






### here i combine two figures together

pdf("../figures/combined_plots.pdf", width = 7.63, height = 48)

plots_combined <- list()
for (i in 1:21) {
  # Create combined plot with plots[i] and plots_ncER[i]
  combined_plot <- plot_grid(plots[[i]], plots_ncER[[i]], ncol = 1, align = 'v')
  plots_combined[[i]] <- combined_plot
}

# Arrange all combined plots in one column
arranged_plots_combined <- plot_grid(plotlist = plots_combined, ncol = 1)
print(arranged_plots_combined)
dev.off()



################################################################################
###prepare motif location talbes with targeted gene, based on category_8,9,10###
################################################################################

# options(scipen = 999)
# prepare the category 8, 9, 10 for centogene, GEL and Rachel, compared to NSC_ESC_common
#  I also added the qval information in order to keep the significant motifs

NSC_whole_table <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>% 
  .[,c(1,8:10)]; colnames(NSC_whole_table) <- c("location", "NSC_TargetGene", "NSC_OMIM_Phenotypes", "NSC_HPO")

ESC_whole_table <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>%
  .[,c(1,8:10)]; colnames(ESC_whole_table) <- c("location", "ESC_TargetGene", "ESC_OMIM_Phenotypes", "ESC_HPO")

# enhancer_motif excel
excel_file <- "../data/patterns/Category_8/enhancer_motif_category8.xlsx"
# excel_file <- "../data/patterns/Category_9/enhancer_motif_category9.xlsx"
# excel_file <- "../data/patterns/Category_10_high/enhancer_motif_category10_high.xlsx"
# excel_file <- "../data/patterns/Category_10_common/enhancer_motif_category10_common.xlsx"

sheets <- excel_sheets(excel_file)
enhancer_motif <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(enhancer_motif) <- sheets

# tomtom_result excel
excel_file <- "../data/patterns/Category_8/tomtom_results_category8.xlsx"
# excel_file <- "../data/patterns/Category_9/tomtom_results_category9.xlsx"
# excel_file <- "../data/patterns/Category_10_high/tomtom_results_category10_high.xlsx"
# excel_file <- "../data/patterns/Category_10_common/tomtom_results_category10_common.xlsx"

sheets <- excel_sheets(excel_file)
tomtom_result <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(tomtom_result) <- sheets

# motif -> GeneName table
motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")

# convert MotifName to GeneName
tomtom_GeneName <- tomtom_result %>% lapply( function(x) tidyr::separate(x, col = pattern, into = c("pattern_type", "pattern_name"), sep = "\\.")) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match0" = "Motif"))) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match1" = "Motif"))) %>%
  lapply( function(x) dplyr::left_join(x, motif2gene, by=c("match2" = "Motif"))) %>%
  lapply( function(x) dplyr::select(x, c("pattern_type", "pattern_name", "num_seqlets", 
                                         "GeneName.x", "GeneName.y", "GeneName", 
                                         "pval0", "pval1", "pval2",
                                         "qval0", "qval1", "qval2"))) %>%
  lapply( function(x) setNames(x, c("pattern_type", "pattern_name", "num_seqlets", 
                                    "motif_GeneName0", "motif_GeneName1", "motif_GeneName2",
                                    "pval0", "pval1", "pval2",
                                    "qval0", "qval1", "qval2")))

# left_join enhancer_motif with Motif's GeneName  
enhancer_motifGeneName <- Map(function(x, y) dplyr::left_join(x, y, by=c("pattern_type", "pattern_name")), enhancer_motif, tomtom_GeneName)

enhancer_motifGeneName <- enhancer_motifGeneName %>% lapply( function(x) dplyr::left_join(x, NSC_whole_table, "location", relationship = "many-to-many" )) %>%
  lapply( function(x) dplyr::left_join(x, ESC_whole_table, "location", relationship = "many-to-many" )) %>%
  lapply( function(x) x[,c(1,15:20,2:14)]) %>%
  lapply( function(x) unique(x)) 

# save as excel
library(openxlsx)
# excel_file <- "../data/patterns/Category_8/enhancer_motif_GeneName_category8.xlsx"
# excel_file <- "../data/patterns/Category_9/enhancer_motif_GeneName_category9.xlsx"
# excel_file <- "../data/patterns/Category_10_high/enhancer_motif_GeneName_category10_high.xlsx"
excel_file <- "../data/patterns/Category_10_common/enhancer_motif_GeneName_category10_common.xlsx"

wb <- createWorkbook()

for (i in seq_along(enhancer_motifGeneName)) {
  sheet_name <- names(enhancer_motifGeneName)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, enhancer_motifGeneName[[i]])
}

saveWorkbook(wb, excel_file, overwrite = T)


# I will further filter those motif_location only with positive_pattern and qval < 0.05 at least in motif0, motif1 or motif2
# also I will combine category8,9,10 at one table
category8 <- read_excel("../data/patterns/Category_8/enhancer_motif_GeneName_category8.xlsx", sheet = "NSC_category8",
                        col_types = c(rep("text",10), "numeric", rep("text",3), rep("numeric", 6))) %>% 
  mutate(Category = "NSC_category8")
category9 <- read_excel("../data/patterns/Category_9/enhancer_motif_GeneName_category9.xlsx", sheet = "NSC_category9",
                        col_types = c(rep("text",10), "numeric", rep("text",3), rep("numeric", 6))) %>% 
  mutate(Category = "NSC_category9")
  
category10_high <- read_excel("../data/patterns/Category_10_high/enhancer_motif_GeneName_category10_high.xlsx", sheet = "NSC_category10",
                         col_types = c(rep("text",10), "numeric", rep("text",3), rep("numeric", 6))) %>% 
  mutate(Category = "NSC_category10")

category10_common <- read_excel("../data/patterns/Category_10_common/enhancer_motif_GeneName_category10_common.xlsx", sheet = "NSC_category10",
                              col_types = c(rep("text",10), "numeric", rep("text",3), rep("numeric", 6))) %>% 
  mutate(Category = "NSC_category10")

Categories <- do.call("rbind", list(category8, category9, category10_high, category10_common))
write.csv(Categories, "../data/patterns/For_WGS_screen/category_8_9_10_motif.csv", row.names = F)

new_categories <- Categories %>% subset(pattern_type == "pos_patterns") %>%
  subset(qval0 < 0.05 | qval1< 0.05 | qval2 < 0.05) %>%
  .[,c(1:4,8:21)] %>%
  unique() %>%
  arrange(location)
write.csv(new_categories, "../data/patterns/For_WGS_screen/category_8_9_10_positive_motif_qval0.05.csv", row.names = F)


new_categories <- Categories %>% subset(pattern_type == "pos_patterns") %>%
  subset(qval0 < 0.05 | qval1< 0.05 | qval2 < 0.05) %>%
  .[,10] %>%
  unique() %>%
  separate(col = motif_location, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  arrange(Chr, as.numeric(Start))
write.table(new_categories, "../data/patterns/For_WGS_screen/category_8_9_10_positive_motif_qval0.05.txt", quote = F, row.names = F, sep = "\t")



#######################################################
###prepare motif location talbes with targeted gene ###
###        based on whole enhancer dataset          ###
###       NSC category_1,2,3,4,5,6,7,8,9,10         ###
#######################################################

# options(scipen = 999)
# prepare the category 1,2,3,4,5,6,7,8, 9, 10 for centogene, GEL and Rachel, compared to NSC_ESC_common
#  I also added the qval information in order to keep the significant motifs

# load NSC enhancer information
NSC_whole_table <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Original_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>% 
  .[,c(1,8:10)]; colnames(NSC_whole_table) <- c("location", "NSC_TargetGene", "NSC_OMIM_Phenotypes", "NSC_HPO")

# load enhancer_motif excels
file_path <- "../data/patterns/Categories/motif_data/"
file_pattern <- "enhancer_motif_category"
file_list <- list.files(path = file_path, pattern = file_pattern, full.names = TRUE)

read_1st_tab <- function(file_path) {
  sheet_name <- excel_sheets(file_path)[1]  # Select the first sheet
  read_excel(file_path, sheet = sheet_name)
} # Function to read only the first sheet of each Excel file and use sheet names as list names

enhancer_motif <- lapply(file_list, read_1st_tab) # Use lapply to read only the first sheet of each Excel file
names(enhancer_motif) <- sapply(file_list, function(file_path) excel_sheets(file_path)[1]) # Set list names as the first sheet names

# load tomtom_result excel
file_path <- "../data/patterns/Categories/motif_data/"
file_pattern <- "tomtom_results_category"
file_list <- list.files(path = file_path, pattern = file_pattern, full.names = TRUE)

tomtom_result <- lapply(file_list, read_1st_tab) # Use lapply to read only the first sheet of each Excel file
names(tomtom_result) <- sapply(file_list, function(file_path) excel_sheets(file_path)[1]) # Set list names as the first sheet names

motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")
tomtom_GeneName <- lapply(tomtom_result, function(x) {
  tidyr::separate(x, col = pattern, into = c("pattern_type", "pattern_name"), sep = "\\.") %>%
    dplyr::left_join(motif2gene, by=c("match0" = "Motif")) %>%
    dplyr::left_join(motif2gene, by=c("match1" = "Motif")) %>%
    dplyr::left_join(motif2gene, by=c("match2" = "Motif")) %>%
    dplyr::select(c("pattern_type", "pattern_name", "num_seqlets", 
                       "GeneName.x", "GeneName.y", "GeneName", 
                       "pval0", "pval1", "pval2",
                       "qval0", "qval1", "qval2")) %>%
    setNames(c("pattern_type", "pattern_name", "num_seqlets", 
                  "motif_GeneName0", "motif_GeneName1", "motif_GeneName2",
                  "pval0", "pval1", "pval2",
                  "qval0", "qval1", "qval2")) 
}) # convert MotifName to GeneName

enhancer_motifGeneName <- Map(function(x, y) dplyr::left_join(x, y, by=c("pattern_type", "pattern_name")), enhancer_motif, tomtom_GeneName)

enhancer_motifGeneName <- Map(
  function(x, list_name) {
    left_join(x, NSC_whole_table, "location", relationship = "many-to-many") %>%
      select(c(1, 15:17, 2:14)) %>%
      mutate(Category = list_name) %>%
      unique()
  },
  enhancer_motifGeneName,
  list_name = names(enhancer_motifGeneName)
)

# combine categories together
Categories <- do.call("rbind", enhancer_motifGeneName)

excel_file <- "../data/patterns/Categories/enhancer_motif_GeneName_categories.xlsx"
write.xlsx(Categories, file = excel_file, append = FALSE, overwrite = T)

##################################################################
###GWAS epilepsy paper, find snps in high-confidence motif site###
##################################################################
options(scipen = 999)
library(data.table)
library(readxl)
library(dplyr)
library(tidyverse)

# convert epilepsy_snp from hg19 to hg38
epilepsy_snp <- read_excel("../data/patterns/GWAS_epilepsy/epilepsy_snp_hg19.xlsx")
epilepsy_snp$chr <- sub("^", "chr", epilepsy_snp$chr)

epilepsy_coordinates <- epilepsy_snp %>% .[,c(4,5)] %>% mutate(Start = pos - 1) %>% .[,c(1,3,2)]

colnames(epilepsy_coordinates) <- c("Chr","Start","End")
epilepsy_coordinates <- epilepsy_coordinates %>% 
  arrange(Chr,as.numeric(Start)) %>%
  unique()

write.table(epilepsy_coordinates, "../data/patterns/GWAS_epilepsy/epilepsy_snp_coordinates.hg19.txt", 
            sep = "\t", row.names = F, col.names = F, quote = F)

# intersect
# hc_motif <- fread("../data/patterns/For_WGS_screen/category_8_9_10_motif.csv") %>%
#   subset(pattern_type == "pos_patterns") %>%
#   mutate(location_tmp = motif_location) %>%
#   separate(col = location_tmp, into = c("chr", "start_tmp"), sep = ":") %>%
#   separate(col = start_tmp, into = c("start", "end"), sep = "-")

hc_motif <- fread("../data/patterns/For_WGS_screen/category_8_9_10_positive_motif_qval0.05.csv") %>%
  mutate(location_tmp = motif_location) %>%
  separate(col = location_tmp, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start)) %>%
  mutate(end = as.numeric(end)) 

epilepsy_snp_hg38 <- read.table("../data/patterns/GWAS_epilepsy/epilepsy_snp_coordinates.hg38.txt", sep = "\t") %>%
  .[,c(1,3,4)] %>%
  separate(col = V4, into = c("chr_hg19", "pos_tmp"), sep = ":") %>%
  separate(col = pos_tmp, into = c("pos_hg19", "End"), sep = "-") %>% 
  .[,c(1,2,4)]
colnames(epilepsy_snp_hg38) <- c("chr", "pos_hg38", "pos_hg19")


# intersect motif_site with snp
result <- epilepsy_snp_hg38 %>%
  left_join(hc_motif, by = "chr", relationship = "many-to-many") %>%
  mutate(intersect = ifelse(pos_hg38 >= start & pos_hg38 <= end, motif_location, NA)) %>%
  select(-start, -end)
result <- result[!is.na(result$intersect), ]

#########################################################
#####GWAS from different studies (External reference)#####
#########################################################
### 1. intersected with enhancers
# load enhancers
enhancers <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Version1/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx") %>%
  select(enhancers, `10_Categories`) %>%
  mutate(location_tmp = enhancers) %>%
  separate(col = location_tmp, into = c("chr_enhancer", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start_enhancer", "end_enhancer"), sep = "-") %>%
  mutate(start_enhancer = as.numeric(start_enhancer)) %>%
  mutate(end_enhancer = as.numeric(end_enhancer)) %>%
  unique()

# load snps
excel_file <- "../External_reference/SNPs.xlsx"
sheets <- excel_sheets(excel_file)
sheets <- sheets[-1] # remove the summary tab
snps <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(snps) <- sheets

# intersect
result <- lapply(snps, function(x) {
  left_join(x, enhancers, by = c("chr_hg38" = "chr_enhancer"), relationship = "many-to-many") %>%
    mutate(intersect = ifelse(pos_hg38 >= start_enhancer & pos_hg38 <= end_enhancer, enhancers, NA)) %>%
    select(-start_enhancer, -end_enhancer) %>%
    filter(!is.na(intersect)) %>% 
    unique()
})
  

# save the intersected snps
library(openxlsx)
excel_file <- "../External_reference/SNPs_enhancers.xlsx"
wb <- createWorkbook()

for (i in seq_along(result)) {
  sheet_name <- names(result)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, result[[i]])
}

saveWorkbook(wb, excel_file, overwrite = T)

### 2. intersected with motifs
# load motifs
pos_motif <- read_excel("../data/patterns/Categories/enhancer_motif_GeneName_categories.xlsx") %>%
  # filter(pattern_type %in% "pos_patterns") %>%
  mutate(location_tmp = motif_location) %>%
  separate(col = location_tmp, into = c("chr_motif", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start_motif", "end_motif"), sep = "-") %>%
  mutate(start_motif = as.numeric(start_motif)) %>%
  mutate(end_motif = as.numeric(end_motif)) 

# intersect
result <- lapply(snps, function(x) {
  
  left_join(x, pos_motif, by = c("chr_hg38" = "chr_motif"), relationship = "many-to-many") %>%
    mutate(intersect = ifelse(pos_hg38 >= start_motif & pos_hg38 <= end_motif, location, NA)) %>%
    select(-start_motif, -end_motif) %>%
    filter(!is.na(intersect)) %>% 
    unique() 
})

# save the intersected snps
library(openxlsx)
excel_file <- "../External_reference/SNPs_motifs.xlsx"
wb <- createWorkbook()


for (i in seq_along(result)) {
  sheet_name <- names(result)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, result[[i]])
}

saveWorkbook(wb, excel_file, overwrite = T)




########################
#####STARR-seq FACS#####
########################

enhancers_activity <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", sheet = 1) %>%
  .[,1:2] %>% unique()

selected_enhancers <- read_excel("../data/patterns/For_WGS_screen/Selected validated enhancers.hg38.xlsx", sheet = "Sheet2")

selected_enhancers <- selected_enhancers %>% left_join(enhancers_activity, by=c('enhancer_location'= 'enhancers'), relationship = "many-to-many") %>%
  arrange(desc(`log2(avg.NSC.RPP+1)`))

####
#### 1. Selected 15 enhancers FACS
####
FACS <- read_excel("../data/Experiments/new/Selected_15_enhancers.xlsx", sheet = 2) %>% 
  gather(key = "Replicate", value = "percent", -c("Group", "Condition"))%>% 
  na.omit()
FACS$Group <- gsub(" .*", "", FACS$Group)

summary_stats <- FACS %>%
  group_by(Group, Condition) %>%
  summarise(mean = mean(percent),
            sd = sd(percent))

FACS <- left_join(FACS, summary_stats, by = c("Group", "Condition"))

FACS$Group <- factor(FACS$Group,
                     levels = c("empty_vector", "OAT", "LAMB2", "ADAR",
                                "ACTB", "TKT","PAFAH1B1", "IRF2BPL", "ZBTB11", "DNMT3A",
                                "NAT8L", "NAA20", "ASH1L", "ATP6V1A", "CIC", "TRIO"))

ggplot(FACS, aes(x = Group, y = mean, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width=0.8) +
  labs(x = "",
       y = "FACS intensity %") +
  scale_fill_discrete() +  # Use discrete color scale
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 100)) +
  geom_errorbar(aes(x = Group, ymin = mean - sd, ymax = mean + sd),
                stat = "identity", position = "dodge", width = 0.3) +
  scale_fill_manual(values = c("#00ba38", "#619cff"))

empty_vector <- FACS$percent[FACS$Group == "empty_vector"]; ACTB <- FACS$percent[FACS$Group == "NAA20"]
(t.test(empty_vector, ACTB, var.equal = T, alternative = "two.sided", paired = F))$p.value



#f16a70

####
#### 2. Selected 15 enhancers SDM FACS
####
FACS <- read_excel("../data/Experiments/new/Selected_15_enhancers_SDM.xlsx", sheet = 2) %>% 
  gather(key = "Replicate", value = "percent", -c("Group", "Condition")) %>% 
  na.omit()
FACS$Group <- gsub(" .*", "", FACS$Group)

summary_stats <- FACS %>%
  group_by(Group, Condition) %>%
  summarise(mean = mean(percent),
            sd = sd(percent))

FACS <- left_join(FACS, summary_stats, by = c("Group", "Condition"))

FACS$Group <- factor(FACS$Group,
                     levels = c("empty_vector", "OAT", "LAMB2", "ADAR",
                                "ACTB", "TKT","PAFAH1B1", "IRF2BPL", "ZBTB11", "DNMT3A",
                                "NAT8L", "NAA20", "ASH1L", "ATP6V1A", "CIC", "TRIO"))
FACS$Condition <- factor(FACS$Condition,
                     levels = c("vector", "WT", "MUT"))

ggplot(FACS, aes(x = Group, y = mean, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width=0.8) +
  labs(x = "",
       y = "FACS intensity %") +
  scale_fill_discrete() +  # Use discrete color scale
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 100)) +
  geom_errorbar(aes(x = Group, ymin = mean - sd, ymax = mean + sd),
                stat = "identity", position = position_dodge(0.8), width = 0.3) +
  scale_fill_manual(values = c("#00ba38", "#619cff", "#f16a70"))

####
#### 3. Core motifs SDM FACS
####
FACS <- read_excel("../data/Experiments/new/FACS-NSC-(SDM)-STARR-seq.xlsx", sheet = 2) %>% 
  gather(key = "Replicate", value = "percent", -c("Group", "Condition")) %>% 
  na.omit()
FACS$Group <- gsub(" .*", "", FACS$Group)

summary_stats <- FACS %>%
  group_by(Group, Condition) %>%
  summarise(mean = mean(percent),
            sd = sd(percent))

FACS <- left_join(FACS, summary_stats, by = c("Group", "Condition"))

# GEL data
GEL <- c("empty_vector", "MN1", "KPTN", "RAB7A", "GRIA4", "ZEB2", "ZEB2_27bpDel", "ZEB2_SNV")
FACS <- subset(FACS, (Group %in% GEL))
FACS$Group <- factor(FACS$Group, levels = GEL)

# # nucleotide convertion: core sites
# core_converstion <- c("empty_vector", "OAT", "ACTB", "PAFAH1B1", "ASH1L")
# FACS <- subset(FACS, (Group %in% core_converstion))
# FACS$Group <- factor(FACS$Group,levels = core_converstion)


FACS$Condition <- factor(FACS$Condition,
                         levels = c("vector", "WT", "MUT"))
# Assuming FACS is the data frame with mean and standard deviation columns
ggplot(FACS, aes(x = Group, y = mean, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width=0.8) +
  labs(x = "",
       y = "FACS intensity %") +
  scale_fill_discrete() +  # Use discrete color scale
  theme_classic() +
  theme(panel.grid=element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 100)) +
  geom_errorbar(aes(x = Group, ymin = mean - sd, ymax = mean + sd),
                stat = "identity", position = position_dodge(0.8), width = 0.3) +
  scale_fill_manual(values = c("#00ba38", "#619cff", "#f16a70"))


# pair-wise
GRIA4_WT <- FACS$percent[FACS$Group == "KPTN WT"]; GRIA4_MUT <- FACS$percent[FACS$Group == "KPTN MUT"]
(t.test(GRIA4_WT, GRIA4_MUT, var.equal = T, alternative = "two.sided", paired = T))$p.value


#########################################
###External paper SNPs, SNP percentile###
#########################################

### 1. percentile for positive and negative SNPs
SNP_percentile <- read.csv("../External_reference/SNP_percentile/SNP_percentile.all_rank.csv") %>% 
  subset(group %in% c("diff", "daSNVs", "emVars")) %>%
  na.omit() 

SNP_percentile <- read.csv("../External_reference/SNP_percentile/SNP_percentile.all_rank.csv") %>% 
  subset(group %in% c("diff", "daSNVs", "emVars")) %>%
  na.omit() %>%
  subset(Percentile > 60)


SNP_percentile$positive <- factor(SNP_percentile$positive, levels = c("yes", "no"))
SNP_percentile$group <- factor(SNP_percentile$group, levels = c("diff", "daSNVs", "emVars"))

ggplot(SNP_percentile, aes(x=group, y=Percentile, fill=positive)) + 
  geom_boxplot(outlier.size=0.2) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("normalize_cb percentile %") +
  scale_fill_manual(values = c("#f8766d", "#619cff"))

diff_yes <- subset(SNP_percentile, (group %in% "diff") & (positive %in% "yes")); diff_no <- subset(SNP_percentile, (group %in% "diff") & (positive %in% "no"))
daSNVs_yes <- subset(SNP_percentile, (group %in% "daSNVs") & (positive %in% "yes")); daSNVs_no <- subset(SNP_percentile, (group %in% "daSNVs") & (positive %in% "no"))
emVars_yes <- subset(SNP_percentile, (group %in% "emVars") & (positive %in% "yes")); emVars_no <- subset(SNP_percentile, (group %in% "emVars") & (positive %in% "no"))
yes <- subset(SNP_percentile, positive %in% "yes"); no <- subset(SNP_percentile, positive %in% "no")

pairwise_results <- c((t.test(diff_yes$Percentile, diff_no$Percentile, var.equal = F, alternative = "greater", paired = F))$p.value,
                      (t.test(daSNVs_yes$Percentile, daSNVs_no$Percentile, var.equal = F, alternative = "greater", paired = F))$p.value,
                      (t.test(emVars_yes$Percentile, emVars_no$Percentile, var.equal = F, alternative = "greater", paired = F))$p.value)
pairwise_results

### 2. percentile for positive and negative SNPs, separate analysis!!!
### correlate the cb percentile with the MPRA FDR ranking

## (1.1) diffSNP: correlation between yes and no, rank by MPRA_FDR
# read the original SNP activity
diffSNP <- read_excel("../External_reference/data_source/5173 schizophrenia gwas snps 439 showing differential activity wt vs disease allele Systematic investigation of allelic regulatory activity of schizophrenia-associated common variants.xlsx") %>%
  mutate(Rank = rank(-MPRA_FDR), FDR_Percentile = ((Rank) / (n())) * 100)
diffSNP$CHR <- sub("^", "chr", diffSNP$CHR)

# read the SNP cb_percentile
SNP_percentile <- read.csv("../External_reference/SNP_percentile/SNP_percentile.each_rank.csv") %>% 
  subset(group %in% "diff") %>%
  separate(col = hg19, into = c("chr_hg19", "pos_tmp"), sep = ":") %>%
  separate(col = pos_tmp, into = c("pos_hg19", "pos_hg19_tmp"), sep = "-") %>%
  mutate(pos_hg19 = as.numeric(pos_hg19)) %>%
  left_join(diffSNP[,c(1,2,5,14)], by=c('chr_hg19'= 'CHR', 'pos_hg19' = 'BP'), relationship = "many-to-many")


SNP_percentile$positive <- factor(SNP_percentile$positive, levels = c("yes", "no"))
ggscatter(SNP_percentile, x = "FDR_Percentile", y = "Percentile",
          color = "positive",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "darkgray", fill = "lightgray")) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = c("yes" = "red", "no" = "blue")) +
  xlab("MPRA_FDR percentile %") +
  ylab("cb percentile %")

## (1.2) diffSNP: correlation between yes and no, subset MPRA_FDR < 0.1 and rank by logFC
# read the original SNP activity
diffSNP <- read_excel("../External_reference/data_source/5173 schizophrenia gwas snps 439 showing differential activity wt vs disease allele Systematic investigation of allelic regulatory activity of schizophrenia-associated common variants.xlsx") %>%
  subset(MPRA_FDR < 0.1) %>%
  mutate(Rank = rank(abs(MPRA_logFC)), FC_Percentile = ((Rank) / (n())) * 100)
diffSNP$CHR <- sub("^", "chr", diffSNP$CHR)

# read the SNP cb_percentile
SNP_percentile <- read.csv("../External_reference/SNP_percentile/SNP_percentile.each_rank.csv") %>% 
  subset((group %in% "diff") & (positive %in% "yes")) %>%
  separate(col = hg19, into = c("chr_hg19", "pos_tmp"), sep = ":") %>%
  separate(col = pos_tmp, into = c("pos_hg19", "pos_hg19_tmp"), sep = "-") %>%
  mutate(pos_hg19 = as.numeric(pos_hg19)) %>%
  left_join(diffSNP[,c(1,2,5,14)], by=c('chr_hg19'= 'CHR', 'pos_hg19' = 'BP'), relationship = "many-to-many")

ggscatter(SNP_percentile, x = "FC_Percentile", y = "Percentile",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "darkgray", fill = "lightgray")) +
  stat_cor(method = "pearson") +
  xlab("MPRA_FC percentile %") +
  ylab("cb percentile %")


## (2.1) daSNVs: mim fdr
# read the original SNP activity
. 

SNP_percentile$positive <- factor(SNP_percentile$positive, levels = c("yes", "no"))
ggscatter(SNP_percentile, x = "FDR_Percentile", y = "Percentile",
          color = "positive",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "darkgray", fill = "lightgray")) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = c("yes" = "red", "no" = "blue")) +
  xlab("MPRA_FDR percentile %") +
  ylab("cb percentile %")


## (2.2) daSNVs: mim fdr
# read the original SNP activity
excel_file <- "../External_reference/data_source/892 daSNVs whole snps with fdr and logFC.xlsx"
sheets <- excel_sheets(excel_file)
sheets <- sheets[-c(1,15,16,17)] # remove the summary tab
daSNVs <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(daSNVs) <- sheets

daSNVs <- lapply(daSNVs, function(x) {
  subset(x, fdr < 0.05) %>%
  separate(col = rowname, into = c("chr_hg19", "pos_hg19"), sep = "_") %>%
    mutate(pos_hg19 = as.numeric(pos_hg19)) %>%
    mutate(logFC = as.numeric(logFC)) %>%
    select("chr_hg19", "pos_hg19", "logFC")
})

merged_daSNVs <- Reduce(function(x, y) merge(x, y, by = c("chr_hg19", "pos_hg19"), all = TRUE), daSNVs)
names(merged_daSNVs)[3:ncol(merged_daSNVs)] <- names(daSNVs)
max_abs_FC <- apply(merged_daSNVs[, 3:ncol(merged_daSNVs)], 1, function(row) max(abs(row), na.rm = TRUE))
merged_daSNVs$max_abs_FC <- max_abs_FC
merged_daSNVs <- mutate(merged_daSNVs, Rank = rank(max_abs_FC), FC_Percentile = ((Rank) / (n())) * 100)


# read the SNP cb_percentile
SNP_percentile <- read.csv("../External_reference/SNP_percentile/SNP_percentile.each_rank.csv") %>% 
  subset(group %in% "daSNVs") %>%
  separate(col = hg19, into = c("chr_hg19", "pos_tmp"), sep = ":") %>%
  separate(col = pos_tmp, into = c("pos_hg19", "pos_hg19_tmp"), sep = "-") %>%
  mutate(pos_hg19 = as.numeric(pos_hg19))  %>%
  left_join(merged_daSNVs, by = c("chr_hg19", "pos_hg19"))

SNP_percentile$positive <- factor(SNP_percentile$positive, levels = c("yes", "no"))
ggscatter(SNP_percentile, x = "FC_Percentile", y = "Percentile",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "darkgray", fill = "lightgray")) +
  stat_cor(method = "pearson") +
  xlab("MPRA_FC percentile %") +
  ylab("cb percentile %")



## (3.1) daSNVs cell types
# read the original SNP activity
excel_file <- "../External_reference/data_source/892 daSNVs whole snps with fdr and logFC.xlsx"
sheets <- excel_sheets(excel_file)
sheets <- sheets[-c(1,15,16,17)] # remove the summary tab
daSNVs <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(daSNVs) <- sheets

daSNVs <- lapply(daSNVs, function(x) {
  mutate(x, Rank = rank(-fdr), FDR_Percentile = ((Rank) / (n())) * 100) %>%
    separate(col = rowname, into = c("chr_hg19", "pos_hg19"), sep = "_") %>%
    mutate(pos_hg19 = as.numeric(pos_hg19)) 
})

# read the SNP cb_percentile
SNP_percentile <- read.csv("../External_reference/SNP_percentile/SNP_percentile.csv") %>% 
  subset(group %in% "daSNVs") %>%
  separate(col = hg19, into = c("chr_hg19", "pos_tmp"), sep = ":") %>%
  separate(col = pos_tmp, into = c("pos_hg19", "pos_hg19_tmp"), sep = "-") %>%
  mutate(pos_hg19 = as.numeric(pos_hg19))  

daSNVs <- lapply(daSNVs, function(df) {
  left_join(df, SNP_percentile, by = c("chr_hg19", "pos_hg19"))
})


# Loop through each data frame in daSNVs and create ggscatter plots
for (i in seq_along(daSNVs)) {
  plot_data <- daSNVs[[i]]
  
  ggscatter(plot_data, x = "FDR_Percentile", y = "Percentile",
            color = "positive",
            add = "reg.line",
            conf.int = TRUE,
            add.params = list(color = "darkgray", fill = "lightgray")) +
    stat_cor(method = "pearson") +
    scale_color_manual(values = c("yes" = "red", "no" = "blue")) +
    ggtitle(paste("Plot for", names(daSNVs)[i]))
  
  # Add any other customization or saving the plot if needed
  # For example, you can save each plot using ggsave function
  ggsave(paste0("plot_", names(daSNVs)[i], ".png"), plot = last_plot())
}

###########################
###daSNVs ratio boxplots###
###########################

# finally only keep rs200483 and rs301806 !!!!
# we only keep SNPs in their MPRA and our plasmids

#### our enhancers
enhancers <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/Version1/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx") %>%
  dplyr::select(enhancers, `10_Categories`) %>%
  mutate(location_tmp = enhancers) %>%
  separate(col = location_tmp, into = c("chr_enhancer", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start_enhancer", "end_enhancer"), sep = "-") %>%
  mutate(start_enhancer = as.numeric(start_enhancer)) %>%
  mutate(end_enhancer = as.numeric(end_enhancer)) %>%
  unique() 
####

#### LD SNPs
LD_snps <- fread("../External_reference/SNP_percentile/examples_guo/GWAS_LD/rs200483.tsv") %>% .[,c(1,2,9)] %>%
  separate(col = "Location", into = c("chr_hg38", "pos_tmp"), sep = ":") %>%
  separate(col = "pos_tmp", into = c("pos_hg38", "ex"), sep = "-") %>%
  .[,-4] %>%
  mutate(pos_hg38 = as.numeric(pos_hg38))
LD_snps$chr_hg38 <- sub("^", "chr", LD_snps$chr_hg38 )
####

raw_count <- read_excel("../External_reference/SNP_percentile/examples_guo/raw_count_boxplots.xlsx", sheet = "counts") %>% 
  .[,c(1:29,48)] %>%
  gather(key = "tissue", value = "rna", -c(hg19, plasmid, chr_hg38, pos_hg38)) %>%
  mutate(rna_dna = rna / plasmid) %>% subset(!(rna_dna %in% c("Inf", "NaN", 0)))
raw_count$tissue <- sub("_.*", "", raw_count$tissue)

sub_count <- raw_count %>%
  left_join(LD_snps, by = c("chr_hg38", "pos_hg38")) %>% # keep LD SNPs
  na.omit() %>% # keep LD SNPs
  left_join(enhancers[,c(1,3:5)], by=c('chr_hg38'= 'chr_enhancer'), relationship = "many-to-many") %>% # keep SNPs in our enhancers
  filter(pos_hg38 <  end_enhancer & pos_hg38 > start_enhancer ) %>% # keep SNPs in our enhancers
  group_by(tissue) %>%
  mutate(zscore = abs((rna_dna - mean(rna_dna)) / sd(rna_dna))) %>%
  filter(zscore < 2.5) %>%
  mutate(condition = ifelse(grepl("Ref", hg19), "Ref", ifelse(grepl("alt", hg19), "Alt", NA))) %>%
  mutate(mean_ref_sign = median(rna_dna[condition == 'Ref'])) %>%
  mutate(ratio = rna_dna / mean_ref_sign) %>%
  ungroup()

# boxplot of MPRA
plot <-sub_count
plot$tissue <- factor(plot$tissue, levels = c("AST", "ES", "N-D2", "N-D4", "N-D10", "A-NPC", "P-NPC", "HEK293"))
plot$condition <- factor(plot$condition, levels = c("Ref", "Alt"))
# plot$`Variation ID` <- factor(plot$`Variation ID`, levels = c("rs301805", "rs301806", "rs301807", "rs301802")) # rs301806
plot$`Variation ID` <- factor(plot$`Variation ID`, levels = c("rs200483", "rs17751184", "rs390764", "rs401754", "rs2747054", "rs200497", "rs67101035")) # rs200483

pdf("../External_reference/SNP_percentile/examples_guo/figures/rs200483_mpra.pdf")
ggplot(plot, aes(x=`Variation ID`, y=ratio, fill=condition)) + 
  geom_boxplot(outlier.size=0.2) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("ratio") +
  scale_fill_manual(values = c("#619cff", "#f8766d")) +
  facet_wrap(~tissue, scales = "free_y") +
  ylim(0, 10)
dev.off()
# normalize_cb 
rm(plot)


####
scores <- read.table("../External_reference/SNP_percentile/examples_guo/scores.txt", sep = "\t", header = T)
plot <- subset(scores, group == "group1")
order <- c("rs301805", "rs301806", "rs301807", "rs301802") # rs301806

# plot <- subset(scores, group == "group2")
# order <- c("rs200483", "rs17751184", "rs390764", "rs401754", "rs2747054", "rs200497", "rs67101035") # rs200483
pdf("../External_reference/SNP_percentile/examples_guo/figures/rs200483_cbDot.pdf", width = 4, height = 8)
ggdotchart(plot, x = "snp", y = "normalize_cb",
           palette = "jco", order=order,
           color = "snp", group = "snp",
           add.params = list(color = "lightgray", size = 1.5),
           dot.size = 15,
           add = "segment", ggtheme = theme_pubr(),
           label = round(plot$normalize_cb, 2), 
           font.label = list(color = "white", size = 15, 
                             vjust = 0.5)) +
  xlab("")
dev.off()
####





# SNPs in enhancers and mpra data
enhancer_snp <- sub_count %>% 
  dplyr::select(c( "pos_hg38", "enhancers", `Variation ID`)) %>%
  unique()
enhancer_ID <- unique(enhancer_snp$enhancers)
pos_ID <- unique(enhancer_snp$pos_hg38)
snp_ID <- unique(enhancer_snp$`Variation ID`)
# LD snps
LD_pos_ID <- LD_snps$pos_hg38
######


###### contribution scores

library(gridExtra)
cb <- fread("../External_reference/SNP_percentile/Whole_contri_percentile.txt")

# one example
plots <- list()
for (i in 1:length(enhancer_ID) ) {
  
  plot <- cb %>% subset(Enhancer %in% enhancer_ID[i])

  plots[[i]] <- ggplot(plot, aes(x = Pos, y = normalize_cb)) +
    geom_line(size = 0.3) +
    geom_point(data = subset(plot, Pos %in% pos_ID), aes(x = Pos, y = normalize_cb), color = "red", size = 0.5) +
    geom_point(data = subset(plot, (Pos %in% LD_pos_ID) & !(Pos %in% pos_ID)), aes(x = Pos, y = normalize_cb), color = "blue", size = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    ggtitle(paste0(enhancer_ID[i], " ", snp_ID[i])) +
    xlab("") +
    theme_bw() +
    theme(panel.grid=element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.5),
          axis.text=element_text(size=8,colour = "black"),
          legend.position="none",
          plot.title = element_text(color="black", size=8)) 
}
pdf("../External_reference/SNP_percentile/examples_guo/figures/rs174561_cb.pdf", width = 7.63, height = 9)

arranged_plots <- do.call(grid.arrange, c(plots[1:8], ncol = 1))

dev.off()
######

###############
###############
###############

#### LD SNPs
LD_snps <- fread("../External_reference/SNP_percentile/examples_guo/GWAS_LD/rs301806.tsv") %>% .[,c(1,2,9)] %>%
  separate(col = "Location", into = c("chr_hg38", "pos_tmp"), sep = ":") %>%
  separate(col = "pos_tmp", into = c("pos_hg38", "ex"), sep = "-") %>%
  .[,-4] %>%
  mutate(pos_hg38 = as.numeric(pos_hg38))
LD_snps$chr_hg38 <- sub("^", "chr", LD_snps$chr_hg38 )
####

snps <- c("chr22.42670965", "chr3.186793242", "chr22.42672124",
          "chr6.27805255", "chr6.27774824",
          "chr1.8482078", "chr1.8481016", "chr1.8504421", "chr1.8484823")

# snps <- c('chr1.8481016', 'chr1.8482078','chr1.8484529','chr1.8484823','chr1.8510577','chr1.8526142')

raw_count <- read_excel("../External_reference/SNP_percentile/examples_guo/raw_count_boxplots.xlsx", sheet = "counts") %>% 
  .[,c(1:29,48)] %>%
  gather(key = "tissue", value = "rna", -c(hg19, plasmid, chr_hg38, pos_hg38)) %>%
  mutate(rna_dna = rna / plasmid) %>% subset(!(rna_dna %in% c("Inf", "NaN", 0)))
raw_count$tissue <- sub("_.*", "", raw_count$tissue)


sub_count <- raw_count %>%
  filter(grepl(paste(snps, collapse="|"), index)) %>%
  group_by(tissue) %>%
  mutate(zscore = abs((rna_dna - mean(rna_dna)) / sd(rna_dna))) %>%
  filter(zscore < 2.5) %>%
  mutate(condition = ifelse(grepl("Ref", index), "Ref", ifelse(grepl("alt", index), "Alt", NA))) %>%
  mutate(mean_ref_sign = median(rna_dna[condition == 'Ref'])) %>%
  mutate(ratio = rna_dna / mean_ref_sign)

sub_count$index <- sub("^([^\\.]+\\.[^\\.]+)\\..*", "\\1",  sub_count$index)
# sub_count$tissue <- sub("_.*", "", sub_count$tissue)

plot <- subset(sub_count, index == "chr1.8484823")
plot$tissue <- factor(plot$tissue, levels = c("AST", "ES", "N-D2", "N-D4", "N-D10", "A-NPC", "P-NPC", "HEK293"))
plot$condition <- factor(plot$condition, levels = c("Ref", "Alt"))
ggplot(plot, aes(x=tissue, y=ratio, fill=condition)) + 
  geom_boxplot(outlier.size=0.2) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("") +
  ylab("ratio") +
  scale_fill_manual(values = c("#619cff", "#f8766d")) 





for (snp in snps) {
  plot <- subset(sub_count, index == snp)
  plot$tissue <- factor(plot$tissue, levels = c("AST", "ES", "N-D2", "N-D4", "N-D10", "A-NPC", "P-NPC", "HEK293"))
  plot$condition <- factor(plot$condition, levels = c("Ref", "Alt"))
  
  p<- ggplot(plot, aes(x = tissue, y = ratio, fill = condition)) + 
    geom_boxplot(outlier.size = 0.2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") +
    ylab("ratio") +
    scale_fill_manual(values = c("#619cff", "#f8766d")) +
    ggtitle(paste("SNP: ", snp))  
  
  ggsave(paste0("snp_plot_", snp, ".png"), plot = p)
}



# plot selected SNPs scores
scores <- read.table("../External_reference/SNP_percentile/examples_guo/scores.txt", sep = "\t", header = T)
scores[is.na(scores)] <- 0 # replace NA as 0 for cb
scores <- scores %>% gather(key = "type", value = "value", -c(snp, group)) 

scores$group <- factor(scores$group, levels = c("group1", "group2", "group3"))
scores$type <- factor(scores$type, levels = c("normalize_cb", "ncER", "CADD"))
scores$snp <- factor(scores$snp, levels = c("rs762995", "rs2239612", "rs134882",
                                             "rs200483", "rs34706883",
                                             "rs301806", "rs301807", "rs301805", "rs159963"))

ggdotchart(scores, x = "group", y = "value", 
           group = "snp", color = "snp", palette = "jco",
           add = "segment", position = position_dodge(0.5),
           facet.by = "type")+
  facet_wrap(~type, scales = "free_y") +
  xlab("")

ggdotchart(scores, x = "snp", y = "value", 
           group = "group", color = "snp", palette = "jco",
           add = "segment", position = position_dodge(0.5),
           facet.by = "type")+
  facet_grid(type ~ group, scales = "free_y") +
  xlab("")


################################################################
###Nijmegen denovo data, find snps in the whole enhancer list###
################################################################
denovo_Nij <- read_excel("../data/patterns/For_WGS_screen/Nijmegen_dnv/20231006_DNMs unsolved cohort_overlap.xlsx") %>%
  .[,c(2,3,5,6,1)] %>% unique()
colnames(denovo_Nij) <- c("chr", "pos", "ref", "alt", "trio_ID")

enhancers <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx")

calculation_DNVs <- function(category_value) {
  category <- enhancers %>%
    # subset(`5_Categories` %in% category_value) %>%
    subset(`10_Categories` %in% category_value) %>%
    dplyr::select(enhancers) %>%
    unique() %>%
    separate(col = enhancers, into = c("chr", "start_tmp"), sep = ":") %>%
    separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.numeric(start)) %>%
    mutate(end = as.numeric(end)) %>%
    unique()
  
  denovo_Nij_enhancer <- denovo_Nij %>%
    left_join(category, by = "chr", relationship = "many-to-many") %>%
    mutate(hit = ifelse(pos >= start & pos <= end, pos, NA)) %>%
    filter(!(is.na(hit))) %>%
    unique()
  
  return(nrow(denovo_Nij_enhancer))
}

calculation_DNVs("Category_10")


#####################################
#####Burden test for gnomAD, GEL#####
#####################################

# Read data
burden_dat <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/burden_test.xlsx", sheet = 1)

# Get unique groups
unique_groups <- unique(burden_dat$group)

# Define subsets for each group
subsets <- lapply(unique_groups, function(group_name) {
  subset(burden_dat, group == group_name)
})


# Initialize an empty list to store results
results_list <- list()

# Process each subset
for (i in seq_along(subsets)) {
  data <- subsets[[i]]
  
  # Check if the subset has only one row
  if (nrow(data) <= 1) {
    next  # Skip to the next subset if there's only one row
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

# Combine all results into a single data frame
combined_results <- do.call(rbind, results_list)

#### plot gnomAD (three diagrams, NSC, ESC and diff_activity)
diagram <- c("codingExon", "random", "threeUTR", "fiveUTR", "DAE", "nDAE",
             "NSC_category_1", "NSC_category_2", "NSC_category_3", 
             "NSC_category_4", "NSC_category_5", "NSC_top10") # NSC in genomAD

diagram <- c("codingExon", "random", "threeUTR", "fiveUTR", "DAE", "nDAE",
             "ESC_category_1", "ESC_category_2", "ESC_category_3", 
             "ESC_category_4", "ESC_category_5", "ESC_top10") # ESC in genomAD

diagram <- c("codingExon", "random", "threeUTR", "fiveUTR", "DAE", "nDAE",
             "common_high", "ESC_high", "NSC_high") # diff_activity in genomAD

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

#### GEL plot figure in groups
#### plot gnomAD (three diagrams, NSC, ESC and diff_activity)
diagram <- c("NSC_category_1", "NSC_category_2", "NSC_category_3", 
             "NSC_category_4", "NSC_category_5", "NSC_top10") # NSC in GEL

diagram <- c("ESC_category_1", "ESC_category_2", "ESC_category_3", 
             "ESC_category_4", "ESC_category_5", "ESC_top10") # ESC in GEL

diagram <- c("common_high", "ESC_high", "NSC_high") # diff_activity in GEL

diagram <- c("common_high", "ESC_high", "NSC_high",
             "ESC_category_1", "ESC_category_2", "ESC_category_3", 
             "ESC_category_4", "ESC_category_5", "ESC_top10",
             "NSC_category_1", "NSC_category_2", "NSC_category_3", 
             "NSC_category_4", "NSC_category_5", "NSC_top10") # all in genomAD

plot_dat <- subset(combined_results, group %in% c("NDD", "control") & category %in% diagram) 
gnomAD_OE <- plot_dat[plot_dat$group == "gnomAD", "odds_ratio"]

plot_dat$category <- factor(plot_dat$category, levels = diagram)
ggplot(plot_dat, aes(y = category, x = odds_ratio, group = group, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +  
  geom_errorbarh(position = position_dodge(width = 0.5), aes(xmin = lower_ci, xmax = upper_ci), height = 0) +
  xlim(c(0.9, 1.4)) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
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
  geom_text(aes(label = sprintf("%.6f", odds_ratio)), hjust = 0.6, vjust = -1.2, size = 2.5) +
  scale_color_manual(values=c("#0E4c92", "#ff6666"))






#### plot GEL cb percentile
diagram <- c("NSC_category_1.60th","NSC_category_2.60th","NSC_category_3.60th","NSC_category_4.60th","NSC_category_5.60th","NSC_top10.60th",
             "NSC_category_1.70th","NSC_category_2.70th","NSC_category_3.70th","NSC_category_4.70th","NSC_category_5.70th","NSC_top10.70th",
             "NSC_category_1.80th","NSC_category_2.80th","NSC_category_3.80th","NSC_category_4.80th","NSC_category_5.80th","NSC_top10.80th",
             "NSC_category_1.90th","NSC_category_2.90th","NSC_category_3.90th","NSC_category_4.90th","NSC_category_5.90th","NSC_top10.90th",
             "NSC_category_1.95th","NSC_category_2.95th","NSC_category_3.95th","NSC_category_4.95th","NSC_category_5.95th","NSC_top10.95th",
             "NSC_category_1.99th","NSC_category_2.99th","NSC_category_3.99th","NSC_category_4.99th","NSC_category_5.99th","NSC_top10.99th") # NSC in genomAD

plot_dat <- subset(combined_results, group == "cb" & category %in% diagram) 
gnomAD_OE <- plot_dat[plot_dat$group == "cb", "odds_ratio"]

plot_dat$category <- factor(plot_dat$category, levels = diagram)
ggplot(plot_dat, aes(y = category, x = odds_ratio)) +
  geom_point(size = 2, color = "#0E4c92") +  
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0, color = "#0E4c92") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 0.5, alpha = 0.5) +
  # scale_y_continuous(name = "", breaks=1:nrow(plot_dat), labels = plot_dat$category, trans = "reverse") +
  xlim(c(0.34,2.1)) + 
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


########################################################
##### enhancers distribution across all chromosomes#####
########################################################

enhancers <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx")  %>%
  .[,1] %>%
  separate(col = enhancers, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>%
  unique() %>%
  subset(chr %in% c(paste0("chr", 1:22), "chrX")) %>%
  mutate(Width = as.numeric(end) - as.numeric(start)) %>%
  group_by(chr) %>%
  count(chr)

# the ratio of (Size_enhancer / Size_chromosome)
chr_distribution <-  read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/burden_test.xlsx", sheet = 4) %>%
  mutate(ratio_Nenhancer=as.numeric(n_enhancer)/as.numeric(chr_size)) %>%
  mutate(ratio_Senhancer=as.numeric(enhancer_size)/as.numeric(chr_size))

chr_distribution$chr <- factor(chr_distribution$chr, levels = c(paste0("chr", 1:22), "chrX"))
ggplot(chr_distribution, aes(x = chr, y = ratio_Senhancer)) +
  geom_bar(stat = "identity") +
  labs(x = "Chromosome", y = "the ratio of (Size_enhancer / Size_chromosome)", title = "enhancer distribution") +
  theme_classic() 

# the number of enhancers
ggplot(chr_distribution, aes(x = chr, y = n_enhancer)) +
  geom_bar(stat = "identity") +
  labs(x = "Chromosome", y = "the number of enhancers", title = "enhancer distribution") +
  theme_classic() 

########################################################
#####YY1 degron data with YY1 high-confidence motif#####
########################################################
YY1_motif <- fread("../data/patterns/For_WGS_screen/category_8_9_10_positive_motif_qval0.05.csv") %>%
  .[,-c(12:17)] %>%
  mutate(location_tmp = motif_location) %>%
  separate(col = location_tmp, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start_motif", "end_motif"), sep = "-") %>% 
  subset(motif_GeneName0 %in% c("ZFP42", "YY2", "Yy1")) %>% 
  filter(!is.na(NSC_OMIM_Phenotypes)) %>% unique()

# NSC YY1 degron enhancer data. We have 0h, 8h_neg, 8h_pos, 24h_neg, 24_pos
excel_file <- "../data/NSC_YY1_RPP/Category_RPP.xlsx"
sheets <- excel_sheets(excel_file)
deg_YY1 <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(deg_YY1) <- sheets

deg_YY1 <- reduce(deg_YY1, left_join, by = "enhancers") %>%
  .[,c(1:7,13,19,25,31)] %>%
  setNames(c("enhancers", 
             "log2(avg_0h+1)", "log2(avg_8h.Neg+1)", "log2(avg_8h.Pos+1)", "log2(avg_24h.Neg+1)", "log2(avg_24h.Pos+1)", 
             "Category_0h", "Category_8h.Neg", "Category_8h.Pos", "Category_24h.Neg", "Category_24.Pos")) %>%
  mutate(enhancers_tmp = enhancers) %>%
  separate(col = enhancers_tmp, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start_deg", "end_deg"), sep = "-")

# intersect motif_site with snp
result <- YY1_motif %>%
  left_join(deg_YY1, by = "chr", relationship = "many-to-many") %>%
  mutate(intersect = ifelse((start_motif >= start_deg & start_motif <= end_deg) |
                              (end_motif >= start_deg & end_motif <= end_deg), 
                            motif_location, NA)) %>% 
  select(-end_deg, -start_deg, -start_motif, -end_motif, -chr) %>%
  filter(!(is.na(intersect))) %>%
  unique()


write.csv(result, "../data/NSC_YY1_RPP/YY1_motif_degron.csv", row.names = F)

lost_motif <- YY1_motif %>% subset(!(motif_location %in% result$motif_location))
write.csv(lost_motif, "../data/NSC_YY1_RPP/YY1_motif_degron_ls.csv", row.names = F)

#################################################
#####word cloud of target gene GO enrichment#####
#################################################
library(wordcloud)
library(openxlsx)

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
    # colors <- brewer.pal(9, "Reds")
    colors <- brewer.pal(9, "Blues")
    pdf(paste0("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Figures/go_enrichment/wordcloud_ESC_", category_names[i], ".pdf"))
    wordcloud(words = go_terms$Term, 
              freq = go_terms[[2]], min.freq = 1, max.words = 200,
              random.order = FALSE, rot.per = 0.35, scale = c(2, 0.01),
              colors = colors)
    dev.off()
  }
}

# Define file path, sheet numbers, and category names
file_path <- "/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 6 - GO BP NSCs and ESCs.xlsx"
# sheet_nums <- 1:6
sheet_nums <- 7:12
category_names <- c("Category_1", "Category_2", "Category_3", "Category_4", "Category_5", "Top10")

# Call the function to create word clouds
create_wordcloud(file_path, sheet_nums, category_names)

#####
##### diff_activity: common_high, NSC_high and ESC_high
#####
go_terms <- read.xlsx("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 7 - GO BP common and cell type specific.xlsx",
                      sheet = 3)
go_terms$Term <- sub("\\s*\\(.*", "", go_terms$Term)
go_terms$log10pvalue <- -log10(go_terms$`P-value`)
go_terms <- go_terms[, c(1, 10)] %>% 
  unique() %>% 
  arrange(desc(log10pvalue)) %>%
  subset(!(Term %in% "regulation of chemokine"))
go_terms$Term <- gsub("negative regulation", "neg. reg.", go_terms$Term)
go_terms$Term <- gsub("positive regulation", "pos. reg.", go_terms$Term)
go_terms$Term <- gsub("regulation", "reg.", go_terms$Term)
go_terms$Term <- gsub("polymerase", "pol", go_terms$Term)
colnames(go_terms) <- c("Term", "log10pvalue")


# colors <- brewer.pal(9, "Reds")
colors <- brewer.pal(9, "Blues")
# colors <- brewer.pal(9, "Greens")

pdf("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Figures/go_enrichment/ESC_high.pdf")
wordcloud(words = go_terms$Term, 
          freq = go_terms[[2]], min.freq = 1, max.words = 200,
          random.order = FALSE, rot.per = 0.35, scale = c(2, 0.01),
          colors = colors)
dev.off()


###############################################
#####categories MotifDiscovery Motif plots#####
###############################################

####
#### prepare category_motif files
####

library(openxlsx)
file_path <- "../data/patterns/Categories/5_categories/motif_data/"
file_names <- list.files(file_path, pattern = "*SC.*\\.txt")
wb <- createWorkbook()

for (file in file_names) {
  
  motif2gene <- read.table("motif_to_gene.txt",sep = "\t"); colnames(motif2gene) <- c("Motif", "GeneName")
  
  full_path <- paste0(file_path, file)  # Construct the full path to the current file
  data <- read.table(full_path, header = F, sep = "\t");
  colnames(data) <- c("pattern",	"num_seqlets",	"match0",	"pval0",	"eval0",	"qval0",
                      "match1",	"pval1",	"eval1",	"qval1",	"match2",	"pval2",	"eval2", "qval2")
  data <- data %>%
    tidyr::separate(col = pattern, into = c("pattern_type", "pattern_name"), sep = "\\.") %>%
    dplyr::left_join(motif2gene, by=c("match0" = "Motif")) %>%
    dplyr::left_join(motif2gene, by=c("match1" = "Motif")) %>%
    dplyr::left_join(motif2gene, by=c("match2" = "Motif")) %>%
    dplyr::select(c("pattern_type", "pattern_name", "num_seqlets", 
                    "GeneName.x", "GeneName.y", "GeneName", 
                    "pval0", "pval1", "pval2",
                    "qval0", "qval1", "qval2")) %>%
    setNames(c("pattern_type", "pattern_name", "num_seqlets", 
               "motif_GeneName0", "motif_GeneName1", "motif_GeneName2",
               "pval0", "pval1", "pval2",
               "qval0", "qval1", "qval2")) 
  
  
  
  sheet_name <- gsub("\\.txt", "", file)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = data)
}

saveWorkbook(wb, file = "../data/patterns/Categories/5_categories/Category_motifs.xlsx", overwrite = TRUE)


####
#### add identified motif_location to enhancers and attach with geneName of motifs
####

# Read all tabs from the Excel file into a list
#####!!!!
# I reorder the tabs in the excel !!!!!
#####!!!!

library(purrr)
excel_file <- "../data/patterns/Categories/5_categories/Category_motifs.xlsx"
all_tabs <- excel_sheets(excel_file)
tabs_to_read <- all_tabs[1:5] 
data_list <- map(tabs_to_read, ~ read_excel(excel_file, sheet = .x) %>%
                   mutate(category = .x))
tomtom_motif <- bind_rows(data_list)
rm(excel_file, all_tabs, tabs_to_read, data_list)

# Directory where your Excel files are located
excel_file <- "../data/patterns/Categories/5_categories/motif_data/"
file_list <- list.files(excel_file, pattern = "^enhancer_motif_NSC.*\\.xlsx$", full.names = TRUE)

data_list <- lapply(file_list, function(file_path) {
  read_excel(file_path, sheet = 1) %>%
    mutate(category = excel_sheets(file_path)[1])
})

motif_enhancer <- bind_rows(data_list)
rm(excel_file, file_list, data_list)

# link with enhancer information with Target gene and HPO terms
enhancers <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx")
  
motif_complete <- motif_enhancer %>% 
  left_join(tomtom_motif, by = c("category", "pattern_type", "pattern_name")) %>%
  .[,c(1,5,2:4,6:15)] %>%
  left_join(enhancers,  by = c("location" = "enhancers"), relationship = "many-to-many") %>%
  .[,c(1,16:24,3:15)]

library(openxlsx)
write.xlsx(motif_complete, "../data/patterns/Categories/5_categories/NSC_5Category_motif.TargetGene.xlsx",colNames = T, rowNames = F)


####
#### Plot category motifs (5_Categories)
####

library(purrr)
excel_file <- "../data/patterns/Categories/5_categories/Category_motifs.xlsx"
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

















####### EnhancerAtlas_Barakat
library(data.table)
motif <- fread("../data/patterns/For_WGS_screen/category_8_9_10_positive_motif_qval0.05.csv")

test <- motif %>%
  separate(col = motif_location, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start_motif", "end_motif"), sep = "-") %>% 
  mutate(start_motif = as.numeric(start_motif) + 5) %>%
  mutate(end_motif = as.numeric(end_motif) - 5) %>%  
  unite("tmp",chr:start_motif, sep = ":") %>%
  unite("motif_location", tmp:end_motif, sep = "-")

write.csv(test, "../data/patterns/For_WGS_screen/category_8_9_10_positive_motif_qval0.05.20bp.csv", row.names = F)

options(scipen = 999)
FB_hg38 <- fread("../data/patterns/For_WGS_screen/FB_enhancers/FB_hg38.txt") %>%
  separate(col = V4, into = c("chr", "start_tmp"), sep = ":") %>%
  separate(col = start_tmp, into = c("start", "end"), sep = "-") %>% 
  mutate(start = as.numeric(start) - 1)


FB_hg19 <- fread("../data/patterns/For_WGS_screen/FB_enhancers/FB_hg19.Gene.bed", colClasses = c("character", "numeric", "character", rep("character", 4))) %>%
  left_join(FB_hg38, by = c("chr", "start", "end"), relationship = "many-to-many") %>%
  filter(!(is.na(V5))) %>%
  .[,c(8:10,4:7)] %>%
  setnames(c("chr", "start", "end","EnhancerType", "GeneID", "GeneName", "GeneDefine"))  %>%
  unite("tmp",chr:start, sep = ":") %>%
  unite("enhancer_location", tmp:end, sep = "-")
write.table(FB_hg19, "../data/patterns/For_WGS_screen/FB_enhancers/FB_hg38.Gene.txt", row.names = F, quote = F, sep = "\t")


# GEL rare variants for genomAD 3.1 AF
RD_phenotype_20bp <- read.csv("../data/patterns/For_WGS_screen/RD_phenotype/RD_phenotype.20bp.csv")
  
variant <- RD_phenotype_20bp %>%  
  mutate(ID = ".") %>%
  mutate(ID1 = ".") %>%
  mutate(ID2 = ".") %>%
  mutate(ID3 = ".") %>%
  .[,c(1,2,19,3,4,20:22)] %>% unique()

write.table(variant, "../data/patterns/For_WGS_screen/RD_phenotype/variants.txt", sep="\t",row.names = F, col.names = F,quote = F)

variants_VEP <- read.table("../data/patterns/For_WGS_screen/RD_phenotype/variants_VEP.txt", sep="\t") %>% 
  .[,c(2,3,27)] %>% unique() %>% 
  separate(col = V2, into = c("chr", "pos_tmp"), sep = ":") %>%
  separate(col = pos_tmp, into = c("pos", "end"), sep = "-") %>% 
  setnames(c("chr", "pos", "end", "alt", "genomAD_AF"))
variants_VEP$chr <- sub("^", "chr", variants_VEP$chr )
variants_VEP$pos <- as.numeric(variants_VEP$pos)
variants_VEP$end <- as.numeric(variants_VEP$end)

variant_gnomAD <- variant[,-c(3,6:8)] %>% left_join(variants_VEP, by = c("chr", "pos"), relationship = "many-to-many") %>%
  mutate(ref_num = apply(.[3],2,nchar)) %>%
  mutate(alt_num.x = apply(.[4],2,nchar)) %>%
  mutate(alt_num.y = apply(.[6],2,nchar)) %>%
  mutate(type = case_when(
    (ref_num == 1 & alt_num.x == 1 & alt.x == alt.y) ~ "substitution",
    (ref_num > 1 & alt_num.x == 1 & alt.y =="-" & ref_num == end -pos +1) ~ "deletion",
    (ref_num == 1 & alt_num.x > 1 & alt_num.x == alt_num.y+1) ~ "insertion")
    ) %>%
  dplyr::filter(!is.na(type)) %>%
  .[,c(1:4,7)] %>%
  unique() %>%
  setnames(c("chr","pos","ref", "alt", "genomAD_AF"))

RD_phenotype_20bp_gnomAD <- variant_gnomAD %>% 
  left_join(RD_phenotype_20bp, by = c("chr", "pos", "ref", "alt"), relationship = "many-to-many")
write.csv(RD_phenotype_20bp_gnomAD, "../data/patterns/For_WGS_screen/RD_phenotype/RD_phenotype.20bp.gnomAD.csv", row.names = F)  



# load snps
excel_file <- "../External_reference/data_source/41588_2023_1533_MOESM8_ESM.xlsx"
sheets <- excel_sheets(excel_file)
sheets <- sheets[-c(14:16)] # remove the summary tab
snps <- lapply(sheets, function(x) read_excel(excel_file, sheet = x))
names(snps) <- sheets

snps <- lapply(snps, function(x) {
  select(x, rowname, fdr, logFC) %>%
    unique()
})

snps <- do.call("rbind", test)

positive <- snps %>% subset(fdr < 0.05)






original <- read_excel("../External_reference/data_source/5173 schizophrenia gwas snps 439 showing differential activity wt vs disease allele Systematic investigation of allelic regulatory activity of schizophrenia-associated common variants.xlsx")
snp <- read_excel("../External_reference/SNPs_comparison.xlsx", sheet = 2) %>%
  separate(col = hg19, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  separate(col = Chr, into = c("pref", "CHR"), sep = 'chr') %>% 
  mutate(Start = as.numeric(Start)) %>%  # Convert 'Start' to numeric
  left_join(original[,c(1,2,9,12)], by=c('CHR'= 'CHR', 'Start' = 'BP'), relationship = "many-to-many") %>%
  mutate(CHR = paste0("chr", CHR)) %>%
  unite("hg19_tmp",CHR:Start, sep = ":") %>%
  unite("hg19", hg19_tmp:End, sep = "-") %>%
  select(-pref) %>%
  arrange(desc(positive), chr_hg38, pos_hg38)
write.table(snp, "../External_reference/439_diff.txt", row.names = F, quote = F, sep = "\t")  


original <- read_excel("../External_reference/data_source/320 snps frVars Functional regulatory variants implicate distinct transcriptional networks in dementia.xlsx", 
                       sheet = 2)
snp <- read_excel("../External_reference/SNPs_comparison.xlsx", sheet = 3) %>%
  separate(col = hg19, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  mutate(Start = as.numeric(Start)) %>%  # Convert 'Start' to numeric
  left_join(original[,c(5,6,11,13)], by=c('Chr'= 'chr', 'Start' = 'pos'), relationship = "many-to-many") %>%
  unite("hg19_tmp",Chr:Start, sep = ":") %>%
  unite("hg19", hg19_tmp:End, sep = "-") %>%
  arrange(desc(positive), chr_hg38, pos_hg38)
write.table(snp, "../External_reference/320_frVars.txt", row.names = F, quote = F, sep = "\t")  


NSC_whole_table <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", 
                              sheet = 2) 
location <- NSC_whole_table %>% 
  subset(`10_Categories` %in% c("Category_10")) %>%
  # subset(!is.na(OMIM_Phenotypes)) %>%
  .[1] %>%
  separate(col = enhancers, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  unique() %>%
  arrange(Chr, as.numeric(Start)) %>%
  mutate(length = as.numeric(End) - as.numeric(Start))

sum(location$length)

write.table(location, "../top10%_OMIM.bed", sep = "\t", quote = F, col.names = F, row.names = F)



NSC_whole_table <- read_excel("/Users/Ruizhi/Work/EMC/Papers/ChIP-STARR-seq NSC/New_version/Supplementary Tables/Supplementary Table 5 - Scaffold, Categories, Targets, HPO.xlsx", 
                              sheet = 1) 
location <- NSC_whole_table %>% 
  subset(`5_Categories` %in% c("Category_4")) %>%
  # subset(`10_Categories` %in% c("Category_10")) %>%
  .[1] %>%
  separate(col = enhancers, into = c("Chr", "Start_tmp"), sep = ":") %>%
  separate(col = Start_tmp, into = c("Start", "End"), sep = "-") %>%
  subset(Chr %in% c(paste0("chr", 1:22), "chrX")) %>%
  unique() %>%
  arrange(Chr, as.numeric(Start)) %>%
  mutate(length = as.numeric(End) - as.numeric(Start) + 1)

sum(location$length)



test <- NSC_whole_table %>%
  # subset(NSC_whole_table$GeneType == "protein-coding") %>%
  select("GeneName") %>% unique()






