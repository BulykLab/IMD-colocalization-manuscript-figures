library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(ggsignif)


### Supp Fig 2A - size and location of OnlyATAC peaks

onlyatac_peaks <- read.table("../data/SuppFig2/Delaneau_OnlyATAC.ATAC.bed")
h3k27ac_peaks <- read.table("../data/SuppFig2/Delaneau_H3K27ac.ATAC.bed")
h3k4me1_peaks <- read.table("../data/SuppFig2/Delaneau_H3K4me1.ATAC.bed")
h3k4me3_peaks <- read.table("../data/SuppFig2/Delaneau_H3K4me3.ATAC.bed")



#ggplot(onlyatac_peaks) + geom_histogram(aes(x=V3-V2), bins = 20) + theme_classic() + xlim(c(0,2000))
#ggplot(h3k27ac_peaks) + geom_histogram(aes(x=V3-V2), bins = 20) + theme_classic() + xlim(c(0,2000))


onlyatac_peaks$Group <- "OnlyATAC"
h3k27ac_peaks$Group <- "H3K27ac"
h3k4me1_peaks$Group <- "H3K4me1"
h3k4me3_peaks$Group <- "H3K4me3"

all_peaks_size <- rbind(onlyatac_peaks, h3k27ac_peaks, h3k4me1_peaks, h3k4me3_peaks)


# wilcox.test(onlyatac_peaks$V3-onlyatac_peaks$V2, h3k27ac_peaks$V3-h3k27ac_peaks$V2) # p < 2.2e-16
# wilcox.test(onlyatac_peaks$V3-onlyatac_peaks$V2, h3k4me1_peaks$V3-h3k4me1_peaks$V2) # p < 2.2e-16
# wilcox.test(onlyatac_peaks$V3-onlyatac_peaks$V2, h3k4me3_peaks$V3-h3k4me3_peaks$V2) # p < 2.2e-16

p_S2A <- all_peaks_size %>% filter(V3-V2 < 10000) %>% ggplot(aes(x=Group, y=V3-V2)) + geom_boxplot(col = "black", fill = "#9F56C8", outlier.alpha = 0.1, outlier.size=0.5) +
  theme_half_open() + background_grid(major = 'y') + geom_hline(yintercept=0) +
  labs(x="Peak set", y="Size of ATAC peaks (bp)") +
  scale_y_continuous(limits = c(0,12000), breaks = c(0, 2500, 5000, 7500, 10000), expand = expansion(mult = c(0, 0))) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        aspect.ratio = 1) +
  geom_signif(
    y_position = c(11600, 11000, 10400), xmin = c(1,2,3), xmax = rep(4,3),
    annotation = rep("*", 3), tip_length = 0.02, textsize = 6, vjust = 0.4
  )


### tss window
histone_peakset_tsswindow <- read.table('~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/IMD_colocalization/draft/data/Fig1/supp/histone_peakset_tsswindow.txt', header = T)


p_S2B <- histone_peakset_tsswindow %>% mutate(tss_window = factor(tss_window, levels=c("10kb", "25kb", "50kb", "100kb", "250kb", "500kb"))) %>%
  ggplot(aes(x=peak_set, group=tss_window)) + theme_classic() + background_grid(major = 'y', minor='y') +
  geom_col( aes(y=num_peaks / Total, alpha=tss_window), color ='black', position = "dodge", fill='#89B2F5') +
  geom_hline(yintercept = 1, color='white') + geom_hline(yintercept = 1, linetype=2) +
  scale_alpha_manual(name="TSS window",
                     labels = c("10kb   (14.5%)", "25kb   (24.7%)", "50kb   (36.7%)", "100kb (51.5%)", "250kb (71.0%)", "500kb (82.8%)"),
                     values = c(0.5,0.6,0.7,0.8,0.9,1)) +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x="Peak set", y="Proportion of ATAC peaks") +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        aspect.ratio = 1)


p_S2 <- p_S2A + plot_spacer() + p_S2B + plot_layout(ncol=3, widths = c(7,1,7))

#p_S1 <- p_S1A / (p_S2B + plot_spacer() + p_S2C + plot_layout(ncol=3, widths = c(7,1,7))) + plot_layout(heights = c(4,5), guides = 'keep')


ggsave('~/Library/CloudStorage/GoogleDrive-rjeong@g.harvard.edu/My Drive/Manuscript/IMD_colocalization/draft/figures/FigureS2AB.pdf', 
       plot = p_S2, 
       width=300, height=120, units="mm")

