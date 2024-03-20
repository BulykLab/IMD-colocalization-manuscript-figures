library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(ggsignif)



## Supplementary Figure 1 - Proportion of h2med by explained by ATAC peaks with various histone marks
peak_histone <- read.table("../data/Fig1/ATAC.mesc.histone_mark.txt", header=T)

p_S1 <- peak_histone %>% mutate(Bin = factor(Bin, levels=c("H3K27ac", "H3K4me1", "H3K4me3", "ELS", "PLS", "OnlyATAC"))) %>%
  ggplot(aes(x=GWAS, group=Bin)) + theme_minimal_hgrid() +
  geom_col( aes(y=h2med_over_h2med_tot, fill = Bin), position = "dodge") +
  geom_errorbar( aes(ymin=h2med_over_h2med_tot-se_h2med_over_h2med_tot, ymax=h2med_over_h2med_tot+se_h2med_over_h2med_tot), position=position_dodge(width = 0.90), width=0.3, colour="black") +
  geom_hline(yintercept = 0) + geom_hline(yintercept = 1, linetype=2) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.position = "right", aspect.ratio = 0.5) +
  background_grid(major = 'y', minor='y') +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#D55E00", "#F0E442", "#CC79A7", "#BBBBBB"), labels=c("H3K27ac (43%)", "H3K4me1 (62%)", "H3K4me3 (15%)", "ELS (30%)", "PLS (14%)", "Only ATAC (35%)")) +
  ylab(expression("Proportion of "*italic(h)["med"]^2)) +
  coord_cartesian(ylim=c(-0.05, 1.8))


ggsave('..//figures/SuppFig1/SuppFigure1.pdf', 
       plot = p_S1, 
       width=250, height=110, units="mm")