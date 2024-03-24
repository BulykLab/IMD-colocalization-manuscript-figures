library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)

source("../code/misc/locuscompare_updated.R")

txdb <- AnnotationDbi::loadDb("../data/supp/txdb_v35_hg38.sqlite")

## Supp Figure 5A - IL6R & CD colocalization

# 1) genome track
chr = "1"
start <- 154236935
end <- 154536935
snp = 'rs9651053'
population = 'EUR'

gr = GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges(start, end))


d=data.frame(x1=154386598, x2=154387307, y1=0.75, y2=1.75)
p_il6r <- ggplot() + theme_classic() +
  geom_rect(data=d, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.6, fill='purple') +
  geom_alignment(
    txdb,
    which = gr,
    cds.rect.h = 0.1,
    color = "black",
    fill = "black",
    label.size = 4,
    arrow.rate = 0,
    length = unit(0.2, "cm"),
    gap.geom = 'segment'
  ) +
  ylim(c(0.75,2.75)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0('chr',chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_point(aes(x = 154386935, y=1), shape=23, size=3, fill="purple") 



# 2) CD GWAS
gwas_file = "../data/SuppFig5/IL6R/IL6R_CD_deLange_2017.sumstat"

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
title = "CD GWAS"
p_gwas <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)



# 3) caQTL
gwas_atac_file = "../data/SuppFig5/IL6R/IL6R.ATAC_chr1_154386598_154387307.Kumasaka.sumstat"

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


# 4) GEUVADIS IL6R eQTL
eqtl_file = "../data/SuppFig5/IL6R/IL6R.LCL_GGGT_meta.sumstat"

merged_eqtl_stat = read.table(eqtl_file, header=T, sep='\t')
merged_eqtl_stat$GEUVADIS_logp = -log10(merged_eqtl_stat$GEUVADIS_p)
merged_eqtl_stat$meta_logp = -log10(merged_eqtl_stat$meta_p)

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
title = "GEUVADIS IL6R eQTL"
p_geuvadis_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)



p_S5A <- p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,7.5)) + background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,14), breaks=c(0,5,10)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,19)) + background_grid(major = 'y') +
  p_il6r + 
  plot_layout(nrow = 4, heights = c(3, 3, 3, 2))

ggsave('../figures/SuppFig5/SuppFig5A.pdf',
       plot = p_S5A, 
       width=130, height=140, units="mm")



## Supp Figure 5C - IL6R & RA no colocalization
# RA GWAS
gwas_file = "../data/SuppFig5/IL6R/IL6R_RA_Ishigaki_2021_EUR.sumstat"

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
p_gwas_ra <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


p_S5C <- p_gwas_ra + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,7.5)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,19)) + background_grid(major = 'y') +
  p_il6r + 
  plot_layout(nrow = 3, heights = c(3, 3, 2))


ggsave('../figures/SuppFig5/SuppFig5C.pdf', 
       plot = p_S5C, 
       width=130, height=90, units="mm")



## Supp Figure 5B - IL12A & PBC colocalization
chr = "3"
start <- 159956116
end <- 160056116
snp = 'rs4679867'
population = 'EUR'

# 1) genome track
gr = GenomicRanges::GRanges(seqnames = "chr3", ranges = IRanges(start, end))

d=data.frame(x1=160005886, x2=160007265, y1=0.75, y2=1.75)
p_il12a <- ggplot() + theme_classic() +
  geom_rect(data=d, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), alpha=0.6, fill='purple') +
  geom_alignment(
    txdb,
    which = gr,
    cds.rect.h = 0.1,
    color = "black",
    fill = "black",
    label.size = 4,
    arrow.rate = 0,
    length = unit(0.2, "cm"),
    gap.geom = 'segment'
  ) +
  ylim(c(0.75,2.75)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0('chr',chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_point(aes(x = 160006116, y=1), shape=23, size=3, fill="purple") 


# 2) PBC GWAS
gwas_file = "../data/SuppFig5/IL12A/IL12A_PBC_Cordell_2021.sumstat"

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
p_gwas <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


# 3) caQTL
gwas_atac_file = "../data/SuppFig5/IL12A/IL12A.ATAC_chr3_160005886_160007265.Kumasaka.sumstat"

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


# 4) GEUVADIS IL12A eQTL
eqtl_file = "../data/SuppFig5/IL12A/IL12A.LCL_GGGT_meta.sumstat"

merged_eqtl_stat = read.table(eqtl_file, header=T, sep='\t')
merged_eqtl_stat$GEUVADIS_logp = -log10(merged_eqtl_stat$GEUVADIS_p)
merged_eqtl_stat$meta_logp = -log10(merged_eqtl_stat$meta_p)

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
title = "GEUVADIS IL12A eQTL"
p_geuvadis_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


p_S5B <- p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,70)) + background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,5), breaks=c(0,2,4)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,20)) + background_grid(major = 'y') +
  p_il12a + 
  plot_layout(nrow = 4, heights = c(3, 3, 3, 2))



ggsave('../figures/SuppFig5/SuppFig5B.pdf', 
       plot = p_S5B, 
       width=130, height=140, units="mm")



## Supp Figure 5D - IL12A & CD no colocalization

## CD GWAS
gwas_file = "../data/SuppFig5/IL12A/IL12A_CD_deLange_2017.sumstat"

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
title = "CD GWAS"
p_gwas_cd <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


p_S5D <- p_gwas_cd + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,7)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,20)) + background_grid(major = 'y') +
  p_il12a + 
  plot_layout(nrow = 3, heights = c(3, 3, 2))

ggsave('../figures/SuppFig5/SuppFig5D.pdf', 
       plot = p_S5D, 
       width=130, height=90, units="mm")
