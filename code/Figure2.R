library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)


## Figure 2A - h2med of ATAC & RNA vs. ATAC vs. RNA
mesc_three <- read.table("../data/Fig2/ATAC_RNA_joint.mesc.results.bytype.070423.txt", header=T)

p_bars_joint <- 
  mesc_three %>%  mutate(QTL = factor(QTL, levels=c("ATAC_RNA", "ATAC", "RNA"))) %>%
  ggplot(aes(x=fct_reorder(GWAS, h2med_prop, .desc=T), group=QTL)) + theme_minimal_hgrid() +
  geom_col( aes(y=h2med_prop, fill= QTL), position = "dodge") +
  geom_errorbar( aes(ymin=h2med_prop-se_h2med_prop, ymax=h2med_prop+se_h2med_prop),position=position_dodge(width = 0.90), width=0.3, colour="black") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_blank(), axis.text.y = element_text(size=14),
        legend.title = element_blank(), legend.text = element_text(size=14),
        legend.text.align = 0, legend.key.size = unit(0.6, 'cm'),
        aspect.ratio = 0.8) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  background_grid(major = 'y', minor='y') +
  scale_fill_manual(labels = c(expression("caQTL "*textstyle(union())*" eQTL"),"caQTL", "eQTL"), values = c("#9F56C8","#89B2F5", "#FF6666")) +
  coord_cartesian(ylim=c(-0.02, 0.6)) +
  ylab(expression(italic(h)["med"]^{2}*" / "*italic(h)["SNP"]^2))
  
p_tile_new <- ggplot(mesc_three, aes(x=fct_reorder(GWAS, h2med_prop, .desc=T), group=QTL)) + theme_void() +
  geom_tile(aes(y=-0.5, fill=Type), color = 'white', height=0.1, width=0.9, show.legend = T) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size=12, angle=45, vjust = 1, hjust=1), axis.text.y = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size=14), legend.key.size = unit(0.6, 'cm'),
        aspect.ratio = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c("#8F4A29", "#BBBBBB"))

p_2A <- p_bars_joint+p_tile_new + plot_layout(nrow=2,guides = 'collect')

ggsave('../figures/Fig2/Figure2A.pdf', 
       plot = p_2A, 
       width=150, height=90, units="mm")


## Figure 2C - Colocalization
loci_count_geuvadis <- read.table('../data/Fig2/caqtl_geuvadis_2.gwas_pw.loci_count.062323.txt', header = T)


p_2C <- loci_count_geuvadis %>% 
  mutate(IMD = factor(IMD, levels=c("PBC", "MS", "RA", "JIA", "SLE", "CD","UC", "IBD", "ATD", "VIT", "AST", "ALL"))) %>%
  ggplot(aes(x=IMD, y=n_loci/Total, fill= Type)) + geom_col() + 
  geom_text(aes(label = n_loci), position = position_stack(vjust = 0.5), color="white", size=4) + 
  theme_classic() +   background_grid(major = 'y') +
  scale_y_continuous(limits = c(0,0.52), expand = expansion(mult = c(0, .02))) +
  scale_fill_manual(labels = c("caQTL & eQTL","Just caQTL", "Just eQTL"), values = c("#9F56C8","#89B2F5", "#FF6666")) +
  labs(y="Proportion of loci colocalized") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_blank(), legend.text = element_text(size=14), legend.key.size = unit(0.5, 'cm'),
        aspect.ratio = 0.5)



ggsave('../figures/Fig2/Figure2C.pdf', 
       plot = p_2C, 
       width=180, height=100, units="mm")


### GO enrichment of ATAC only peaks
great_enrichment <- read.table('../data/Fig2/GREAT_caQTL_only_peaks.txt', header = T, sep='\t')

p_2E <- great_enrichment %>% slice(1:13) %>%
  ggplot(aes(x=-log10(pvalue), y=fct_reorder(GO_biological_processes, pvalue, .desc=T))) + geom_col( fill='#89B2F5') +
  theme_classic() +   background_grid(major = 'x') +
  scale_x_continuous(expand = expansion(mult = c(0, .02))) +
  labs(x=expression(-log[10]*"(Enrichment p)"), y="GO Biological Processes") +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        aspect.ratio = 1
  )


ggsave('../figures/Fig2/Figure2E.pdf', 
       plot = p_2E, 
       width=180, height=120, units="mm")