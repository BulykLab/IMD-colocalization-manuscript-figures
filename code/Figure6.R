library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(stringr)
library(AnnotationDbi)
library(ggbio)
library(GenomicRanges)
library(locuscomparer)
library(data.table)

source("../code/supp/locuscompare_updated.R")

## Figure 6A - TNFSF15 colocalization

# 1) genome track
chr = "9"
start = 114660662
end = 114952160
snp = 'rs7848647'
population = 'EUR'

txdb <- AnnotationDbi::loadDb("../data/supp/txdb_v35_hg38.sqlite")

gr = GenomicRanges::GRanges(seqnames = "chr9", ranges = IRanges(start, end))


d=data.frame(x1=114805662, x2=114807160, y1=0.75, y2=1.75)
p_tnfsf15 <- ggplot() + theme_classic() +
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
  ylim(c(0.75,1.75)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x = element_text(size=14), axis.text.x = element_text(size=12)
  ) +ylab("") + xlab(paste0('chr',chr,' (Mb)'))+
  scale_x_continuous(labels=function(x){sprintf('%.2f',x/1e6)}, expand = expansion(mult = c(0, 0)), limit=c(start,end)) +
  geom_point(aes(x = 114806766, y=1), shape=23, size=3, fill="purple") 



# 2) IBD GWAS
gwas_file = "../data/Fig6/TNFSF15/TNFSF15_IBD_deLange_2017.sumstat"

gwas_stat = read.table(gwas_file, header=T, sep='\t')
gwas_stat$logp = -log10(gwas_stat$P)

gwas_stat = gwas_stat[, c('ID', 'logp')]
colnames(gwas_stat)<-c("rsid","logp")

#gwas_stat = read_metal(gwas_file, marker_col = 'rsid', pval_col = 'pval')
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
gwas_atac_file = "../data/Fig6/TNFSF15/TNFSF15.ATAC_chr9_114805662_114807160.Kumasaka.sumstat"

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




# 4) BLUEPRINT monocyte eQTL
eqtl_file = "../data/Fig6/TNFSF15/TNFSF15.BLUEPRINT_monocyte.sumstat"

eqtl_stat = read.table(eqtl_file, header=T, sep='\t')
eqtl_stat$logp = -log10(eqtl_stat$pvalue)

eqtl_stat = eqtl_stat[, c('variant_id', 'logp')]
colnames(eqtl_stat)<-c("rsid","logp")

eqtl_stat = get_position(eqtl_stat, genome = 'hg38')
ld = retrieve_LD(chr, snp, population)
color = assign_color2(eqtl_stat$rsid, snp, ld)

shape = ifelse(eqtl_stat$rsid == snp, 23, 21)
names(shape) = eqtl_stat$rsid

size = ifelse(eqtl_stat$rsid == snp, 4, 2)
names(size) = eqtl_stat$rsid


eqtl_stat$label = ifelse(eqtl_stat$rsid == snp, snp, '')

metal = eqtl_stat[, c('rsid', 'logp', 'chr', 'pos', 'label')]
title = "Monocyte TNFSF15 eQTL"
p_monocyte_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)

# 5) LCL meta TNFSF15 eQTL
eqtl_file = "../data/Fig6/TNFSF15/TNFSF15.LCL_GGGT_meta.sumstat"

merged_eqtl_stat = read.table(eqtl_file, header=T, sep='\t')
merged_eqtl_stat$meta_logp = -log10(merged_eqtl_stat$meta_p)

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
title = "LCL TNFSF15 eQTL"
p_meta_eqtl <- make_locuszoom2_nolabel(metal,title,chr,color,shape,size,range=c(start, end), ylab_linebreak=FALSE)



p_6A_TNFSF15 <- p_gwas + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,30)) + background_grid(major = 'y') +
  p_atac + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,10), breaks = c(0,5,10)) + background_grid(major = 'y') +
  p_meta_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size=10), axis.title.y = element_text(size=14)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,5.5), breaks = c(0,2,4)) + background_grid(major = 'y') +
  p_monocyte_eqtl + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limit=c(0,55), breaks=c(0,20,40)) + background_grid(major = 'y') +
  p_tnfsf15 + 
  plot_layout(nrow = 5, heights = c(3, 3, 3, 3, 2))

p_4C_TNFSF15


ggsave('../figures/Fig6/Figure6A.pdf', 
       plot = p_6A_TNFSF15, 
       width=130, height=140, units="mm")




### Figure 6B - monocyte vs LCL expression levels
mean_tpm_monocyte <- read.table('../data/Fig6/LCL_v_Mono_CPM.txt', header=T)


p_6B <- mean_tpm_monocyte[,c( 'Gene', 'LCL_mean_TPM', 'monocyte_mean_TPM', 'status')] %>% 
  filter(status != "justLCL") %>%
  distinct(.keep_all = TRUE) %>%
  ggplot() + geom_point(aes(x=LCL_mean_TPM+0.01, y=monocyte_mean_TPM+0.01, color=status), size = 2, alpha=0.7) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_vline(xintercept = 1.01, linetype=2) + 
  scale_x_continuous(trans = 'log10') + scale_y_continuous(trans = 'log10') +
  scale_color_manual(values = c("lightgray", "blue", 'blue'), labels=c("LCL & monocyte", "Monocyte only")) +
  labs(color = "IMD-eQTL colocalization", x="LCL mean TPM+0.01", y="Monocyte mean TPM+0.01") +
  theme_minimal_hgrid() +
  geom_point(x=log10(0.29517400+0.01), y=log10(2.226620+0.01), shape=23, size=3.5, fill="red") +
  geom_text(label = "TNFSF15", x=log10(0.29517400), y=log10(2.226620), size=4, hjust=1.2, vjust=0) +
  theme(axis.title.y = element_text(size=14), axis.title.x = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.text = element_text(size=12), legend.title = element_text(size=12), 
        aspect.ratio = 1)

ggsave('../figures/Fig6/Figure6B.pdf', 
       plot = p_6B, 
       width=160, height=120, units="mm")
