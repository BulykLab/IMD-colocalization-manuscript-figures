library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(ggsignif)


### Figure 1

## Figure 1A - Heritability enrichment of ATAC peaks

# Read ATAC-seq peak LDSC results
ldsc_atac <- read.table("../data/Fig1/ATAC.s-ldsc.results.071623.txt", header=T)

# Identifying significant enrichment
ldsc_significant <- ldsc_atac %>% filter(Enrichment_p < 0.05/16)

p_ldsc <- ldsc_atac %>% 
  mutate(GWAS = factor(GWAS, levels=c("PBC", "MS", "CEL", "RA", "JIA", "SLE", "CD","UC", "IBD", "ATD", "VIT", "AST", "ALL", "SCZ", "T2D", "CAD"))) %>%
  ggplot(aes(x=GWAS)) + theme_minimal_hgrid() +
  geom_col( aes(y=Enrichment, fill = Type)) +
  geom_errorbar( aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=0.3, colour="black") +
  geom_hline(yintercept = 1) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16), 
        axis.text.x = element_text(size=12, angle=45, vjust = 0.6),
        axis.text.y = element_text(size=14),
        legend.title = element_blank(), legend.text = element_text(size=14),
        legend.position = "bottom") +
  coord_cartesian(ylim=c(-0.4, 18)) +
    background_grid(major = 'y', minor='y') +
  scale_fill_manual(values = c("#8F4A29", "#FDD262", "#BBBBBB")) +
  ylab("Heritability enrichment") +
  geom_text(data=ldsc_significant, aes(y=Enrichment+Enrichment_std_error+0.5), label = "*", size=6)

## Figure 1B - Mediated heritability of caQTLs

# Read ATAC-seq peak MESC results
mesc_atac <- read.table("../data/Fig1/ATAC.mesc.results.112022.txt", header=T)

mesc_atac$h2med_prop_p <- (1-pnorm(mesc_atac$h2med_prop/mesc_atac$se_h2med_prop))
mesc_atac$h2med_p <- (1-pnorm(mesc_atac$h2med/mesc_atac$se_h2med))

mesc_significant <- mesc_atac %>% filter(h2med_prop_p < 0.05/16)

p_mesc <- mesc_atac %>% 
  mutate(GWAS = factor(GWAS, levels=c("PBC", "MS", "CEL", "RA", "JIA", "SLE", "CD","UC", "IBD", "ATD", "VIT", "AST", "ALL", "SCZ", "T2D", "CAD"))) %>%
  ggplot(aes(x=GWAS)) + theme_minimal_hgrid() +
  geom_col( aes(y=h2med_prop, fill= Type)) +
  geom_errorbar( aes(ymin=h2med_prop-se_h2med_prop, ymax=h2med_prop+se_h2med_prop), width=0.3, colour="black") +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=16), 
        axis.text.x = element_text(size=12, angle=45, vjust = 0.6),
        axis.text.y = element_text(size=14),
        legend.title = element_blank(), legend.text = element_text(size=14),
        legend.position = "bottom") +
  background_grid(major = 'y', minor='y') +
  scale_fill_manual(values = c("#8F4A29", "#FDD262", "#BBBBBB")) +
  ylab(expression(italic(h)["med"]^{2}*" / "*italic(h)["SNP"]^2)) +
  coord_cartesian(ylim=c(-0.012, 0.6)) +
  geom_text(data=mesc_significant, aes(y=h2med_prop+se_h2med_prop+0.02), label = "*", size=6)

p_1AB <- p_ldsc + p_mesc + plot_layout(ncol = 2)

ggsave('../figures/Fig1/Figure1AB.pdf', 
       plot = p_1AB, 
       device='pdf', 
       width=250, height=100, units="mm")




## Figure 1C - Mediated heritability enrichment by peak strength

peak_strength <- read.table("../data/Fig1/ATAC.mesc.peak_strength.txt", header=T)
peak_strength$h2med_enrichment_qvalue <- qvalue(p = peak_strength$h2med_enrichment_pvalue)$qvalues

p_peak_strength <- peak_strength %>% filter((h2med_enrichment_qvalue < 0.05) & (h2med_enrichment > 0)) %>%
  ggplot(aes(x=Bin, y=GWAS)) + geom_point(aes(size=h2med_enrichment, color=-log10(h2med_enrichment_pvalue))) + theme_classic() +
  scale_size_continuous(limits = c(1,10), breaks=c(2,4,6,8,10)) + 
  scale_x_continuous(expand = expansion(mult = c(0.15, .15)), limits=c(1,5))  +  scale_y_discrete(limits=rev) +
  scale_color_continuous(type = "viridis", limits=c(0.6, 6)) +
  labs(x="Peak strength quintile") + 
  theme(axis.title.x = element_text(size=12), axis.title.y = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position='none')

## Figure 1D - Mediated heritability enrichment by histone marks

peak_histone <- read.table("../data/Fig1/ATAC.mesc.histone_mark.txt", header=T)
peak_histone$h2med_enrichment_qvalue <- qvalue(p = peak_histone$h2med_enrichment_pvalue, pi0=0.16)$qvalues


p_mesc_histone <- peak_histone %>% mutate(Bin = factor(Bin, levels=c("H3K27ac", "H3K4me1", "H3K4me3", "ELS", "PLS", "OnlyATAC"))) %>%
  filter((h2med_enrichment_qvalue < 0.05) & (h2med_enrichment > 0)) %>%
  ggplot(aes(x=Bin, y=GWAS)) + geom_point(aes(size=h2med_enrichment, color=-log10(h2med_enrichment_pvalue))) + theme_classic() +
    scale_x_discrete(labels=c("H3K27ac\n(43%)", "H3K4me1\n(62%)", "H3K4me3\n(15%)", "ELS\n(30%)", "PLS\n(14%)", "OnlyATAC\n(35%)")) + scale_y_discrete(limits=rev) +
  scale_color_continuous(type = "viridis", limits=c(0.6, 6)) +
  scale_size_continuous(limits = c(1,10), breaks=c(2,4,6,8,10)) +
  labs(x="Histone mark set", size=expression(italic(h)["med"]^{2}*" enrichment") , color=expression(-log[10]*"(p)")) + #coord_flip() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size=12, angle=45, vjust = 0.6),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=12), legend.text = element_text(size=12))

p_1CD <- p_peak_strength + p_mesc_histone + plot_layout(ncol = 2, widths = c(1,2.2))

ggsave('../figures/Fig1/Figure1CD.pdf', 
       plot = p_1CD, 
       device='pdf', 
       width=200, height=120, units="mm")