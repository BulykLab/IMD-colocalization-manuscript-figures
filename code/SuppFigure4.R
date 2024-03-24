library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)


## Supp Figure 4 - GO Enrichment of IMD-colocalized eQTLs

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

ggsave('../figures/SuppFig4/SuppFig4.pdf', 
       plot = p_S4, 
       width=200, height=100, units="mm")