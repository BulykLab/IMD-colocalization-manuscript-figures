library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)



## Supp Figure 3 - h2med proportion of ATAC & RNA vs. ATAC vs. RNA
mesc_proportion <- read.table('../data/SuppFig3/ATAC_RNA_joint.mesc.results.proportion.070423.txt', header = T)

p_S3 <- mesc_proportion %>% 
  mutate(IMD = factor(IMD, levels=c("PBC", "MS", "RA", "SLE", "CD","UC", "IBD"))) %>%
  ggplot(aes(x=IMD, y=Proportion, fill= QTL)) + geom_col() + 
  geom_text(aes(label = round(Proportion, 2)), position = position_stack(vjust = 0.5), color="white", size=5) + 
    theme_classic() +   background_grid(major = 'y') +
  scale_y_continuous(limits = c(0,0.52), expand = expansion(mult = c(0, .02))) +
  scale_fill_manual(labels = c(expression("caQTL "*intersect()*" eQTL"),"Just caQTL", "Just eQTL"), values = c("#9F56C8","#89B2F5", "#FF6666")) +
  labs(y=expression(italic(h)["med"]^{2}*" / "*italic(h)["SNP"]^2)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_blank(), legend.text = element_text(size=14), legend.text.align = 0,
        legend.key.size = unit(0.6, 'cm'),
        aspect.ratio = 0.75)


ggsave('../figures/SuppFig3/SuppFig3.pdf', 
       plot = p_S3, 
       width=200, height=115, units="mm")
