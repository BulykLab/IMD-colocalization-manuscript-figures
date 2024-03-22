library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(data.table)

source("../code/supp/locuscompare_updated.R")


#### Figure 2A - h2med of ATAC & RNA vs. ATAC vs. RNA
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


#### Figure 2C - Colocalization
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


#### Figure 2D - ELMO1 colocalization plot

chr = "7"
start <- 37242861 
end <- 37442861
snp = 'rs60600003'
population = 'EUR'

txdb <- AnnotationDbi::loadDb("../data/supp/txdb_v35_hg38.sqlite")

gr = GenomicRanges::GRanges(seqnames = "chr7", ranges = IRanges(start, end))

d=data.frame(x1=37340793, x2=37343675, y1=0.75, y2=1.75)
p_elmo1 <- ggplot() + theme_classic() +
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
  geom_point(aes(x=37342861), y=1, shape=23, size=4, fill="purple")


## RA GWAS
gwas_file = "../data/Fig2/ELMO1/ELMO1_RA_Ishigaki_2021_EUR.sumstat"

gwas_stat = read.table(gwas_file, header=T, sep='\t')
gwas_stat$logp = -log10(gwas_stat$P)

gwas_stat = gwas_stat[, c('ID', 'logp')]
colnames(gwas_stat) <- c("rsid","logp")

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
p_ragwas <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


# PBC GWAS
gwas_file = "../data/Fig2/ELMO1/ELMO1_PBC_Cordell_2021.sumstat"

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
title = "PBC GWAS"
p_pbcgwas <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


## MS GWAS
gwas_file = "../data/Fig2/ELMO1/ELMO1_MS_Patsopoulos_2019.sumstat"

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
title = "MS GWAS"
p_msgwas <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


## caQTL
gwas_atac_file = "../data/Fig2/ELMO1/ELMO1.ATAC_chr7_37340793_37343675.Kumasaka.sumstat"

gwas_atac_stat = read.table(gwas_atac_file, header=T, sep='\t')
gwas_atac_stat$atac_logp = -log10(gwas_atac_stat$pval_nominal)

atac_stat = gwas_atac_stat[, c('variant_id', 'atac_logp')]
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


## GEUVADIS & Meta ELMO1 eQTL
eqtl_file = "../data/Fig2/ELMO1/ELMO1.LCL_GGGT_meta.sumstat"

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
title = "LCL ELMO1 eQTL"
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
title = "Meta-analyzed LCL ELMO1 eQTL"
p_meta_eqtl <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


p_2D <- p_ragwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,7.5), breaks = c(0, 2, 4, 6, 8)) + background_grid(major = 'y') +
  p_pbcgwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,15), breaks = c(0, 5, 10,15)) + background_grid(major = 'y') +
  p_msgwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=10), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,8), breaks = c(0, 2, 4, 6, 8)) + background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,7)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,8), breaks = c(0, 2, 4, 6)) + background_grid(major = 'y') +
  p_elmo1 + 
  plot_layout(nrow = 6, heights = c(3, 3, 3, 3, 3, 2))


ggsave('../figures/Fig2/Figure2D.pdf', 
       plot = p_2D, 
       width=120, height=140, units="mm")




#### Figure 2E - GO enrichment of ATAC only peaks
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