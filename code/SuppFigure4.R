library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)


## Supp Figure 4 - GO Enrichment of IMD-colocalized eQTLs
#loci_count_geuvadis <- read.table('~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/IMD_colocalization/draft/data/Fig2/caqtl_geuvadis_2.gwas_pw.loci_count.062323.txt', header = T)
#go_enrichment <- read.table('~/Downloads/GO_enrichment_nearby.txt', header = T, sep='\t')
go_enrichment <- read.table('../data/SuppFig4/GO_eQTL_enrichment_nearby.txt', header = T, sep='\t')


p_S4 <- go_enrichment %>% slice(1:20) %>%
  ggplot(aes(x=-log10(pvalue), y=fct_reorder(GO_biological_processes, pvalue, .desc=T))) + geom_col(fill='#FF6666') +
  theme_classic() +   background_grid(major = 'x') +
  scale_x_continuous(expand = expansion(mult = c(0, .02))) +
  labs(x=expression(-log[10]*"(p)"), y="GO Biological Processes") +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=10),
        aspect.ratio = 1.2
        )
#scale_fill_manual(labels = c("Enrichment","Depletion"), values = c("#9F56C8","#89B2F5")) +  

ggsave('~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/IMD_colocalization/draft/figures/FigureS4.pdf', 
       plot = p_S4, 
       width=200, height=100, units="mm")

aspect.ratio = 0.6

#+ 
  geom_text(aes(label = n_loci), position = position_stack(vjust = 0.5), color="white", size=4) + 
  theme_classic() +   background_grid(major = 'y') +
  scale_y_continuous(limits = c(0,0.52), expand = expansion(mult = c(0, .02))) +
  scale_fill_manual(labels = c("caQTL & eQTL","Just caQTL", "Just eQTL"), values = c("#9F56C8","#89B2F5", "#FF6666")) +
  labs(y="Proportion of loci colocalized") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_blank(), legend.text = element_text(size=14), legend.key.size = unit(0.5, 'cm'),
        aspect.ratio = 0.6)

