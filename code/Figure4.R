library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(stringr)
library(ggsignif)
library(ggtext)
library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(data.table)

source("../code/supp/locuscompare_updated.R")

## Figure 4A - additional eQTL colocalization through meta-anaylsis
loci_gain_meta <- read.table('../data/Fig4/meta_added_eQTL_caQTL_num.070623.test.txt', header = T)

p_4A <- loci_gain_meta %>%
mutate(Type = factor(Type, levels=c("caQTL_only", "LCLmeta", "Geuvadis"))) %>%
ggplot(aes(x=IMD, y=num_loci, fill=Type, linetype=Type)) + geom_col(color='black', width = 0.8) + 
  theme_classic() +  background_grid(major = 'y') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y="Number of\ncaQTL-colocalized loci") + 
  scale_fill_manual(labels = c("no LCL eQTL", "+ LCL meta-analysis", "Geuvadis"), values = c("white", "#FF7600", '#9F56C8')) +
  scale_linetype_manual(labels = c("no LCL eQTL", "+ LCL meta-analysis", "Geuvadis"), values = c('dashed', 'solid', 'solid')) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.text = element_text(size=12), legend.title = element_blank(),
        aspect.ratio = 0.35, legend.position = "none" )


### Figure 4B - TSS distance violin plot

## First, Geuvadis eQTL + caQTL distance
geuvadis_distance <- read.table('../data/Fig4/geuvadis_eQTL_caQTL.distance.070723.txt', header = T)

## Second, meta-analyzed eQTL + caQTL distance
meta_added_distance <- read.table('../data/Fig4/lcl_meta_added.distance.070723.txt', header = T)

geuvadis_distance$Group <- "Geuvadis"
geuvadis_distance$distance <- ifelse((geuvadis_distance$start_atac < geuvadis_distance$end_rna) & (geuvadis_distance$end_atac > geuvadis_distance$end_rna), 0, geuvadis_distance$distance) 
meta_added_distance$Group <- "Meta-analysis"
meta_added_distance$distance <- ifelse((meta_added_distance$start_atac < meta_added_distance$end_rna) & (meta_added_distance$end_atac > meta_added_distance$end_rna), 0, meta_added_distance$distance) 

all_eqtl_caqtl_distance <- rna_atac_h2_split[,c('peak', 'gene', 'peak_distance')]
colnames(all_eqtl_caqtl_distance) <- c('Peak', 'Gene', 'distance')
all_eqtl_caqtl_distance$Group <- "All"

all_distance <- rbind(all_eqtl_caqtl_distance, geuvadis_distance[,c('Peak', 'Gene', 'distance', 'Group')], meta_added_distance[,c('Peak', 'Gene', 'distance', 'Group')])



p_4B <- ggplot(all_distance, aes(x=Group, y=distance+1)) +
  geom_violin(aes(fill=Group, color=Group), width=0.6, alpha=0.5) +
  geom_boxplot(width=0.1) +
  scale_y_continuous(limits=c(1,9e6), trans='log10', breaks = c(1,1e2,1e4,1e6)) +
  labs(y="Peak-to-TSS distance (bp)") +
  scale_fill_manual(values = c("black", '#9F56C8', "#FF7600")) +
  scale_color_manual(values = c("black", '#9F56C8', "#FF7600")) +
  scale_x_discrete(labels = c("All\n(n = 3457)", "Geuvadis\n(n = 90)", "Added by\nmeta-analysis\n(n = 47)")) +
  theme_classic() + background_grid(major = 'y', minor='y') +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        aspect.ratio = 0.5, legend.position="none") +
  geom_signif(
    y_position = c(6.2, 6.7), xmin = c(1.05, 1.05), xmax = c(1.95,2.95),
    annotation = c("*", "*"), tip_length = .02, textsize = 5, vjust = 0.5
  ) +
  geom_signif(
    y_position = c(6.2), xmin = c(2.05), xmax = c(2.95),
    annotation = c("ns"), tip_length = .02, textsize = 4, vjust = -0.1
  )


#wilcox.test(geuvadis_distance$distance, meta_added_distance$distance, alternative = 'less')  # p = 0.05776
#wilcox.test(all_eqtl_caqtl_distance$distance, geuvadis_distance$distance, alternative = 'less')  # p = 0.002283
#wilcox.test(all_eqtl_caqtl_distance$distance, meta_added_distance$distance, alternative = 'less') # p = 5.894e-06

p_4AB <- p_4A + (plot_spacer() +cowplot::get_legend(p_4A+scale_fill_manual(labels = c("Geuvadis", "+ LCL meta-analysis", "no LCL eQTL"), values = c( '#9F56C8', "#FF7600", "white")) +
                                                 scale_linetype_manual(labels = c("Geuvadis", "+ LCL meta-analysis", "no LCL eQTL"), values = c('solid', 'solid', 'dashed')) +
                                                 theme(legend.position = 'top')) + plot_layout(ncol=2, widths=c(1,1e10))) + plot_layout(nrow=2, heights = c(21,1))  +
  plot_spacer() + p_4B + plot_layout(nrow=4, heights = c(21,1,1,30))


ggsave('../figures/Fig4/Figure4AB.pdf', 
       plot = p_4AB, 
       width=160, height=150, units="mm")



### Figure 4C - POU3F1 colocalization plot

chr = "1"
start <- 38028536
end <- 38238536
snp = 'rs7544568'
population = 'EUR'


## 1) genome track
txdb <- AnnotationDbi::loadDb("../data/supp/txdb_v35_hg38.sqlite")

gr = GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(start, end))

d=data.frame(x1=38173053, x2=38174797, y1=0.75, y2=1.75)
p_pou3f1 <- ggplot() + theme_classic() +
  geom_rect(data=d, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.6, fill='purple') +
  geom_alignment(
    txdb,
    which = gr,
    cds.rect.h = 0.1,
    color = "black",
    fill = "black",
    label.size = 5,
    arrow.rate = 0,
    length = unit(0.2, "cm"),
    gap.geom = 'segment'
  ) +
  ylim(c(0.75,1.75)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0('chr',chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_point(aes(x=38173536), y=1, shape=23, size=4, fill="purple")


## 2) RA GWAS
gwas_file = "../data/Fig4/POU3F1/POU3F1_RA_Ishigaki_2021_EUR.sumstat"

gwas_stat = read.table(gwas_file, header=T, sep='\t')
gwas_stat$logp = -log10(gwas_stat$P)

gwas_stat = gwas_stat[, c('ID', 'logp')]
colnames(gwas_stat)<-c("rsid","logp")

gwas_stat = get_position(gwas_stat, genome = 'hg38')

ld = retrieve_LD(chr, snp, population)
color = assign_color2(gwas_stat$rsid, snp, ld)

shape = ifelse(gwas_stat$rsid == snp, 23, 21)
names(shape) = gwas_stat$rsid

size = ifelse(gwas_stat$rsid == snp, 4, 2)
names(size) = gwas_stat$rsid

gwas_stat$label = ifelse(gwas_stat$rsid == snp, snp, '')

metal = gwas_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
title = "RA GWAS"
p_gwas <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


## 3) caQTL
gwas_atac_file = "../data/Fig4/POU3F1/RA_Ishigaki_2021_EUR_ATAC_chr1_38173053_38174797.gwaspw.txt"

gwas_atac_stat = read.table(gwas_atac_file, header=T, sep='\t')
gwas_atac_stat$atac_logp = -log10(pchisq(gwas_atac_stat$Z_ATAC_chr1_38173053_38174797^2, df = 1, lower.tail = F))

atac_stat = gwas_atac_stat[, c('SNPID', 'atac_logp')]
colnames(atac_stat)<-c("rsid","logp")

atac_stat = get_position(atac_stat, genome = 'hg38')
ld = retrieve_LD(chr, snp, population)
color = assign_color2(atac_stat$rsid, snp, ld)

shape = ifelse(atac_stat$rsid == snp, 23, 21)
names(shape) = atac_stat$rsid

size = ifelse(atac_stat$rsid == snp, 4, 2)
names(size) = atac_stat$rsid

atac_stat$label = ifelse(atac_stat$rsid == snp, snp, '')

metal = atac_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
title = "LCL caQTL"
p_atac <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)




## GEUVADIS & Meta POU3F1 eQTL
eqtl_file = "../data/Fig4/POU3F1/POU3F1.LCL_GGGT_meta.sumstat"

merged_eqtl_stat = read.table(eqtl_file, header=T, sep='\t')
merged_eqtl_stat$GEUVADIS_logp = -log10(merged_eqtl_stat$GEUVADIS_p)
merged_eqtl_stat$meta_logp = -log10(merged_eqtl_stat$meta_p)

## First GEUVADIS
geuvadis_eqtl_stat = merged_eqtl_stat[, c('variant_id', 'GEUVADIS_logp')]
colnames(geuvadis_eqtl_stat)<-c("rsid","logp")

eqtl_stat = get_position(geuvadis_eqtl_stat, genome = 'hg38')

ld = retrieve_LD(chr, snp, population)
color = assign_color2(eqtl_stat$rsid, snp, ld)

shape = ifelse(eqtl_stat$rsid == snp, 23, 21)
names(shape) = eqtl_stat$rsid

size = ifelse(eqtl_stat$rsid == snp, 4, 2)
names(size) = eqtl_stat$rsid


eqtl_stat$label = ifelse(eqtl_stat$rsid == snp, snp, '')

metal = eqtl_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
title = "GEUVADIS POU3F1 eQTL"
p_geuvadis_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)

## Then meta-analyzed
meta_eqtl_stat = merged_eqtl_stat[, c('variant_id', 'meta_logp')]
colnames(meta_eqtl_stat)<-c("rsid","logp")

eqtl_stat = get_position(meta_eqtl_stat, genome = 'hg38')
ld = retrieve_LD(chr, snp, population)
color = assign_color2(eqtl_stat$rsid, snp, ld)

shape = ifelse(eqtl_stat$rsid == snp, 23, 21)
names(shape) = eqtl_stat$rsid

size = ifelse(eqtl_stat$rsid == snp, 4, 2)
names(size) = eqtl_stat$rsid


eqtl_stat$label = ifelse(eqtl_stat$rsid == snp, snp, '')

metal = eqtl_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
title = "Meta-analyzed POU3F1 eQTL"
p_meta_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)



p_4C <- p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,11), breaks = c(0, 5, 10)) + background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,23), breaks = c(0, 10,20)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,14), breaks = c(0, 5,10)) + background_grid(major = 'y') +
  p_meta_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,14), breaks = c(0, 5, 10)) + background_grid(major = 'y') +
  p_pou3f1 + 
  plot_layout(nrow = 5, heights = c(3, 3, 3, 3, 2))


ggsave('../figures/Fig4/Figure4C.pdf', 
       plot = p_4C, 
       width=140, height=140, units="mm")
