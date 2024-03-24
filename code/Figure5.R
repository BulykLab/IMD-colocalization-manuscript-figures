library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(stringr)


### Figure 5A - Loci count with immune cell eQTLs

loci_count_immune_celltype <- read.table('../data/Fig5/ATAC_IMD_coloc_immune_loci_count.071023.txt', header = T)

p_5A <- loci_count_immune_celltype %>% filter(Type != "ATAC_only") %>% filter(num > 0) %>%
  mutate(Type = factor(Type, levels=c("Immune", "LCL", "Geuvadis"))) %>%
  ggplot(aes(x=IMD, y=num/total, fill=Type)) + geom_col() + 
  geom_text(aes(label = num), color = "white", position = position_stack(vjust = 0.5), size=5) + 
  theme_classic() +   background_grid(major = 'y') +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values = c("#1CCBD5", "#E46500", "#9F56C8"), labels=c("+ immune cell eQTL", "+ LCL meta-analysis", "Geuvadis") ) +
  labs(y="Proportion of caQTL loci \n with eQTL colocalization") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.position = "bottom", aspect.ratio = 0.7)



### Figure 5B - Loci count with immune cell eQTLs
eqtl_loci_count_immune <- read.table('../data/Fig5/IMD_eQTL_coloc_count_immune_celltype.3.txt', header = T, sep='\t')


p_5B <- eqtl_loci_count_immune %>% 
  mutate(Celltype = factor(Celltype, levels=c("LCL", "B cell", "CD4 T cell", "CD8 T cell", "NK cell", "Macrophage", "Monocyte", "Neutrophil"))) %>%
  mutate(Stimulated = factor(Stimulated, levels=c("naive", "stimulated", "Meta-analyzed"))) %>%
  ggplot() + 
  geom_point(aes(x=Sample_size, y=MS, color=Celltype, shape=Stimulated, fill=Celltype), alpha= 0.8, size=4) +
  theme_classic() +
  scale_x_continuous(trans="log10") + scale_y_continuous(limits=c(0,NA), expand = expansion(mult = c(0, .05))) +
  labs(x="eQTL sample size", y="Number of MS colocalized loci", color="Cell type", fill="Cell type") +
  scale_color_manual(values=c("black", "darkgray", "#5500cc", "#777fff", "#7BE590", "red", "#ff7f00", "#ffd700")) +
  scale_shape_manual(values=c(16,1,23), labels = c("Naive", "Stimulated", "Meta-analyzed")) +
  scale_fill_manual(values=c("black", "darkgray", "#5500cc", "#777fff", "#7BE590", "red", "#ff7f00", "#ffd700")) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_blank(), legend.text = element_text(size=12),
        aspect.ratio = 1)


p_5AB <- p_5A + plot_spacer() + p_5B +plot_layout(ncol = 3, widths = c(10,0.2,7))


ggsave('../figures/Fig5/Figure5AB.pdf', 
       plot = p_5AB, 
       width=290, height=100, units="mm")