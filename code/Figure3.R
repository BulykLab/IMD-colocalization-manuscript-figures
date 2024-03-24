library(ggplot2)
library(forcats)
library(cowplot)
library(patchwork)
library(dplyr)
library(qvalue)
library(stringr)
library(ggsignif)
library(ggtext)



### Figure 3A - Peak-to-TSS distance and cis-h2

###
# Functions

# Copied from https://www.itcodar.com/r/split-violin-plot-with-ggplot2.html
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             # Original function by Jan Gleixner (@jan-glx)
                             # Adjustments by Wouter van der Bijl (@Axeman)
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

###
# Data and Figure
rna_atac_h2 <- read.table("../data/Fig3/ATAC_RNA.h2.atacpos.txt", header=T)

rna_atac_h2$peak_distance <- ifelse((rna_atac_h2$start_atac < rna_atac_h2$gene_pos) & (rna_atac_h2$end_atac > rna_atac_h2$gene_pos), 0, rna_atac_h2$peak_distance) 

#range((rna_atac_h2 %>% filter(peak_distance_quintile == 1))$peak_distance) # 14 5695
#range((rna_atac_h2 %>% filter(peak_distance_quintile == 2))$peak_distance) # 5728 24286
#range((rna_atac_h2 %>% filter(peak_distance_quintile == 3))$peak_distance) # 24300 53647
#range((rna_atac_h2 %>% filter(peak_distance_quintile == 4))$peak_distance) # 53655 111011
#range((rna_atac_h2 %>% filter(peak_distance_quintile == 5))$peak_distance) # 111330 594675


#### TSS distance quintile

rna_h2 <- rna_atac_h2[,c(1:8,12:17)]
colnames(rna_h2)[6:8] <- c("h2cis", "h2cis_se", "h2cis_p")
rna_h2$type <- "RNA"
rna_h2 <- rna_h2[order(rna_h2$peak_distance_quintile, decreasing = FALSE), ]
rna_h2 <- rna_h2[!duplicated(rna_h2$gene), ]
atac_h2 <- rna_atac_h2[,c(1:5,9:17)]
colnames(atac_h2)[6:8] <- c("h2cis", "h2cis_se", "h2cis_p")
atac_h2$type <- "ATAC"
atac_h2 <- atac_h2[order(atac_h2$peak_distance_quintile, decreasing = FALSE), ]
atac_h2 <- atac_h2[!duplicated(atac_h2$peak), ]

rna_atac_h2_split <- rbind(rna_h2, atac_h2)


p_3A <- 
  ggplot(rna_atac_h2_split, aes(x = factor(peak_distance_quintile), y = h2cis, fill = type)) +
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha=1) +
  theme_classic() + background_grid(major = 'y') + 
  labs(x="Peak-to-TSS distance quintile", y=expression(paste("QTL ", italic('h'), {}^2, {}[cis]))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.05)), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_fill_manual(labels=c("caQTL", "eQTL"), values =c("#89B2F5", "#FF6666")) + #"#0F63EB", "#FF4040")) +
  scale_x_discrete(labels=c("1\n(0 ~  \n  5.7 kb)", "2\n(5.7 ~  \n  24.3 kb)", "3\n(24.3 ~  \n  53.6 kb)", "4\n(53.6 ~  \n  111.3 kb)", "5\n(111.3 ~  \n  594.7 kb)")) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_blank(), legend.text = element_text(size=12),
        legend.position = 'right', aspect.ratio = 0.7) +
  geom_signif(
    y_position = rep(-.05,5), xmin = c(0.9, 1.9,2.9,3.9,4.9), xmax = c(1.1, 2.1,3.1,4.1,5.1),
    annotation = c("***", "***", "**", "***", "***"), tip_length = -.02, textsize = 4, vjust = 2
  ) +
  geom_signif(
    y_position = c(1,1.04,1.08,1.12), xmin = rep(1,4), xmax = c(2,3,4,5),
    annotation = c("*", "*", "**", "**"), tip_length = .01, textsize = 4, color="#FF6666", vjust = 0.7
  ) +
  geom_signif(
    y_position = c(1.16,1.2,1.24,1.28), xmin = rep(1,4), xmax = c(2,3,4,5),
    annotation = rep("ns", 4), tip_length = .01, textsize = 3, color="#89B2F5", vjust = 0.2
  )


p_3A


## Wilcoxon test p values
w_pval_one <- rep(NA, 4)
w_pval_two <- rep(NA, 4)
w_pval_three <- rep(NA, 5)

# eQTL h2cis contrast with 1st quintile
for (i in 1:4) {
  w_pval_one[i] <-  wilcox.test((rna_atac_h2_split %>% filter(type=="RNA" & peak_distance_quintile == i+1))$h2cis, (rna_atac_h2_split %>% filter(type=="RNA" & peak_distance_quintile == 1))$h2cis, alternative='l')$p.value
}

# caQTL h2cis contrast with 1st quintile
for (i in 1:4) {
  w_pval_two[i] <-  wilcox.test((rna_atac_h2_split %>% filter(type=="ATAC" & peak_distance_quintile == i+1))$h2cis, (rna_atac_h2_split %>% filter(type=="ATAC" & peak_distance_quintile == 1))$h2cis, alternative='l')$p.value
}

# caQTL vs eQTL h2cis contrast
for (i in 1:5) {
  w_pval_three[i] <-  wilcox.test((rna_atac_h2 %>% filter(peak_distance_quintile == i))$h2cis_ATAC, (rna_atac_h2 %>% filter(peak_distance_quintile == i))$h2cis_RNA, paired = TRUE, alternative = 'g')$p.value
}


### Figure 3B - Regression coefficients for caQTL vs. eQTL h2cis
rna_atac_h2$peak_distance_adjusted <- rna_atac_h2$peak_distance / 100000
lm_h2cis <- as.data.frame(as.matrix(summary(lm(h2cis_RNA ~ h2cis_ATAC + peak_distance_adjusted, rna_atac_h2))$coefficients[2:3,], ncol=4))
lm_h2cis$var <- rownames(lm_h2cis)

p_3B <- ggplot(lm_h2cis, aes(x=var, y=Estimate)) + geom_point(size = 2, color = "black") +
  scale_x_discrete(limits=rev, labels = c("Peak-to-TSS\n   distance\n   (100 kb)", expression(paste("caQTL ", italic('h'), {}^2, {}[cis])))) + #str_wrap("Peak-to-TSS distance (100 kb)", 10))) +
  scale_y_continuous(limits = c(-0.07, 0.15)) +
  geom_errorbar( aes(ymin=Estimate-`Std. Error`, ymax=Estimate+`Std. Error`), width=0.05, colour="black") +
  theme_classic() +  background_grid(major='x') + coord_flip() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12, color='black'),
        aspect.ratio = 0.8) +
  geom_abline(slope = 0, intercept = 0, linetype=1) +
  labs(y="Regression coefficient (SE)")


### Figure 3C - TSS window MESC
tss_window <- read.table("../data/Fig3/ATAC.mesc.tss_window.txt", header=T)


p_3C <- 
  tss_window %>% mutate(Bin = factor(Bin, levels=c("10kb", "25kb", "50kb", "100kb", "250kb", "500kb", "1000kb"))) %>%
    filter(Bin != "50kb") %>%
  ggplot(aes(x=GWAS, group=Bin)) + theme_classic() +
  geom_col( aes(y=h2med_over_h2med_tot, alpha=Bin), position = "dodge", fill='#89B2F5')+  ##526A93') +
  geom_errorbar( aes(ymin=h2med_over_h2med_tot-se_h2med_over_h2med_tot, ymax=h2med_over_h2med_tot+se_h2med_over_h2med_tot), position=position_dodge(width = 0.90), width=0.3, colour="black") +
  geom_hline(yintercept = 1, linetype=2) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
        legend.title = element_text(size=12), legend.text = element_text(size=11),
        legend.position = "right",
        aspect.ratio = 0.6) +
  background_grid(major = 'y', minor='y') +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), expand = expansion(mult = c(0, 0))) +
  scale_alpha_manual(name = "TSS window",values = c(0.4,0.5,0.6,0.7,0.8,0.9,1),
                     labels = c("10 kb   (14.5%)", "25 kb   (24.7%)", "100 kb (51.5%)", "250 kb (71.0%)", "500 kb (82.8%)", "1 Mb    (91.7%)")
  ) +
  ylab(expression("Proportion of "*italic(h)["med; caQTL"]^2))


p_3AB <- p_3A  | (plot_spacer() / p_3B + plot_layout(nrow=2, heights = c(1,10))) 
p_3ABC <- (p_3AB) / (p_3C | plot_spacer()) + plot_layout(heights = c(7,4))




ggsave('../figures/Fig3/Figure3ABC.pdf', 
       plot = p_3ABC, 
       width=300, height=170, units="mm")
