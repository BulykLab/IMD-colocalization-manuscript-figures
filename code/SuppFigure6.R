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

txdb <- AnnotationDbi::loadDb("../data/misc/txdb_v35_hg38.sqlite")

## Supp Figure 6A - CIITA & IBD colocalization

# 1) genome track
chr = "16"
start <- 10731269
end <- 11011269
snp = 'rs12922863'
population = 'EUR'

gr = GenomicRanges::GRanges(seqnames = "chr16", ranges = IRanges(start, end))

d=data.frame(x1=10866552, x2=10872027, y1=0.75, y2=1.75)
p_ciita <- ggplot() + theme_classic() +
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
  ylim(c(0.75,2.75)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0('chr',chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_point(aes(x = 10871269, y=1), shape=23, size=3, fill="purple")


# 2) IBD GWAS
gwas_file = "../data/SuppFig6/CIITA/CIITA_IBD_deLange_2017.sumstat"

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
title = "IBD GWAS"
p_gwas <- make_locuszoom2(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


# 3) caQTL
gwas_atac_file = "../data/SuppFig6/CIITA/CIITA.ATAC_chr16_10866552_10872027.Kumasaka.sumstat"

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



## GEUVADIS & Meta CIITA eQTL
eqtl_file = "../data/SuppFig6/CIITA/CIITA.LCL_GGGT_meta.sumstat"

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
title = "GEUVADIS CIITA eQTL"
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
title = "Meta-analyzed CIITA eQTL"
p_meta_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


p_S6A <- p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,9), breaks = c(0,4,8)) + background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,9), breaks = c(0,4,8)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,9), breaks = c(0,4,8)) + background_grid(major = 'y') +
  p_meta_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,9), breaks = c(0,4,8)) + background_grid(major = 'y') +
  p_ciita + 
  plot_layout(nrow = 5, heights = c(3, 3, 3, 3, 2))


ggsave('../figures/SuppFig6/SuppFig6A.pdf', 
       plot = p_S6A, 
       width=170, height=140, units="mm")


## Supp Figure 6B - ATG16L1 & CD colocalization


chr = "2"
start <- 233047890
end <- 233447890
snp = 'rs11679791'
population = 'EUR'


# 1) genome track
gr = GenomicRanges::GRanges(seqnames = "chr2", ranges = IRanges(start, end))

d=data.frame(x1=233247701, x2=233248180, y1=0.75, y2=1.75)
p_atg16l1 <- ggplot() + theme_classic() +
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
  geom_point(aes(x = 233247890, y=1), shape=23, size=4, fill="purple")


# 2) CD GWAS
gwas_file = "../data/SuppFig6/ATG16L1/ATG16L1_CD_deLange_2017.sumstat"

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
atac_file = "../data/SuppFig6/ATG16L1/ATG16L1.ATAC_chr2_233247701_233248180.Kumasaka.sumstat"

atac_stat = read.table(atac_file, header=T, sep='\t')
atac_stat$atac_logp = -log10(atac_stat$pval_nominal)

atac_stat = atac_stat[, c('variant_id', 'atac_logp')]
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


# 4) GEUVADIS & Meta ATG16L1 eQTL
eqtl_file = "../data/SuppFig6/ATG16L1/ATG16L1.LCL_GGGT_meta.sumstat"

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
title = "GEUVADIS ATG16L1 eQTL"
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
title = "Meta-analyzed ATG16L1 eQTL"
p_meta_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)


p_S6B <- p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_blank()) + 
  background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,10), breaks=c(0,5,10)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,10), breaks=c(0,5,10)) + background_grid(major = 'y') +
  p_meta_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=12), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,10), breaks=c(0,5,10)) + background_grid(major = 'y') +
  p_atg16l1 + 
  plot_layout(nrow = 5, heights = c(3, 3, 3, 3, 2))


ggsave('../figures/SuppFig6/SuppFig6B.pdf', 
       plot = p_S6B, 
       width=130, height=140, units="mm")


## Supp Figure 6C - CARD9 & CD colocalization
chr = "9"
start <- 136323081
end <- 136483081
snp = 'rs4078099'
population = 'EUR'

# 1) genome track
gr = GenomicRanges::GRanges(seqnames = "chr9", ranges = IRanges(start, end))

d=data.frame(x1=136372421, x2=136374081, y1=0.75, y2=1.75)
p_card9 <- ggplot() + theme_classic() +
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
  geom_point(aes(x = 136373081, y=1), shape=23, size=3, fill="purple") 



# 2) CD GWAS
gwas_file = "../data/SuppFig6/CARD9/CARD9_CD_deLange_2017.sumstat"

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


# 2) caQTL
gwas_atac_file = "../data/SuppFig6/CARD9/CARD9.ATAC_chr9_136372421_136374081.Kumasaka.sumstat"

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


# 3) GEUVADIS & Meta CARD9 eQTL
eqtl_file = "../data/SuppFig6/CARD9/CARD9.LCL_GGGT_meta.sumstat"

merged_eqtl_stat = read.table(eqtl_file, header=T, sep='\t')
merged_eqtl_stat$GEUVADIS_logp = -log10(merged_eqtl_stat$GEUVADIS_p)
merged_eqtl_stat$meta_logp = -log10(merged_eqtl_stat$meta_p)

## First GEUVADIS
geuvadis_eqtl_stat = merged_eqtl_stat[, c('variant_id', 'GEUVADIS_logp')]
colnames(geuvadis_eqtl_stat)<-c("rsid","logp")

eqtl_stat = get_position(geuvadis_eqtl_stat, genome = 'hg38')
#snp = 'rs4077515'
ld = retrieve_LD(chr, snp, population)
color = assign_color2(eqtl_stat$rsid, snp, ld)

shape = ifelse(eqtl_stat$rsid == snp, 23, 21)
names(shape) = eqtl_stat$rsid

size = ifelse(eqtl_stat$rsid == snp, 4, 2)
names(size) = eqtl_stat$rsid


eqtl_stat$label = ifelse(eqtl_stat$rsid == snp, snp, '')

metal = eqtl_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
title = "GEUVADIS CARD9 eQTL"
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
title = "Meta-analyzed CARD9 eQTL"
p_meta_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)



p_S6C <- p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,40)) + background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,7)) + background_grid(major = 'y') +
  p_geuvadis_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,24)) + background_grid(major = 'y') +
  p_meta_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=10), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,110)) + background_grid(major = 'y') +
  p_card9 + 
  plot_layout(nrow = 5, heights = c(3, 3, 3, 3, 2))

ggsave('../figures/SuppFig6/SuppFig6C.pdf', 
       plot = p_S6C, 
       width=130, height=140, units="mm")
