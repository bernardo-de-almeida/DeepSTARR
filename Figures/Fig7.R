
library(GenomicRanges)
library(ggplot2)
theme_set(theme_light() + theme(axis.text = element_text(colour = "black")))
library(cowplot)
library(dplyr)

########
# Fig 7A
########

df_twist <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Twist_mutagenesis_Drosophila_final_table.txt", stringsAsFactors = F)[,c(1:5,7,8,27,31)]
df_twist <- df_twist[complete.cases(df_twist$dev_log2FoldChange) & complete.cases(df_twist$hk_log2FoldChange),]

tmp <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Table_synthetic_enhancers.txt")
tmp$obs <- tmp$dev_log2FoldChange
tmp$pred <- tmp$Predictions_dev

### add wt enhancers - just to keep same axis limits
tmp_wt <- df_twist[df_twist$Enhancer_type %in% c("dev") & df_twist$dev_log2FoldChange>3.149511 & df_twist$dev_log2FoldChange>df_twist$hk_log2FoldChange,]
tmp_wt$obs <- tmp_wt$dev_log2FoldChange

tmp$pred_bins <- cut(tmp$pred, breaks = c(floor(min(tmp$pred)), 2,3,4,5, ceiling(max(tmp$pred))), include.lowest = T)

# names
bins_labels <- sapply(levels(tmp$pred_bins), function(x){
  paste0(x, "\n(n=", table(tmp$pred_bins)[x], ")")
})

gg_boxplot <- ggplot() +
  geom_boxplot(data=tmp_wt, aes("Native", obs),
               fill="orangered") +
  geom_boxplot(data=tmp, aes(pred_bins, obs, alpha=pred_bins),
               fill="grey30") +
  scale_alpha_manual(values=seq(0.3,1,length.out = 6)) +
  guides(col=F, alpha=F) +
  scale_y_continuous(paste0("STARR-seq enhancer activity [log2]"),
                     # limits = lim,
                     breaks = seq(-12,20,2)) +
  scale_x_discrete(paste0("DeepSTARR predicted activity [log2]"), limits=c("Native", levels(tmp$pred_bins)),
                   labels=c("Native", bins_labels)) +
  ggtitle(paste0("Synthetic enhancers (n=", nrow(tmp), ")")) +
  theme_bw() +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.grid = element_blank(), axis.line=element_line(colour="black"),
        axis.text=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"),
        axis.title.x=element_text(margin = margin(5,0,0,0)), axis.title.y=element_text(margin = margin(0,5,0,0)),
        plot.title = element_text(size=18, hjust = 0.5, colour="black"), plot.subtitle = element_text(size=14, hjust = 0.5))

gg_scater <- ggplot(tmp, aes(pred, obs)) +
  geom_point(size=1, col="grey30", aes(alpha=pred_bins)) +
  geom_point(data=tmp_wt, aes(1, obs),
             col=NA) +
  scale_alpha_manual(values=seq(0.4,1,length.out = 6)) +
  guides(col=F, alpha=F) +
  scale_y_continuous(paste0("STARR-seq enhancer activity [log2]"),
                     # limits = lim,
                     breaks = seq(-12,20,2)) +
  scale_x_continuous(paste0("DeepSTARR predicted activity [log2]"),
                     # limits = lim,
                     breaks = seq(-12,20,2)) +
  ggtitle(paste0("Synthetic enhancers (n=", nrow(tmp), ")")) +
  theme_bw() +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.grid = element_blank(), axis.line=element_line(colour="black"),
        axis.text=element_text(size=13, colour="black"), axis.title=element_text(size=16, colour="black"),
        axis.title.x=element_text(margin = margin(5,0,0,0)), axis.title.y=element_text(margin = margin(0,5,0,0)),
        plot.title = element_text(size=18, hjust = 0.5, colour="black"), plot.subtitle = element_text(size=14, hjust = 0.5)) +
  annotate("text",  x=min(tmp$pred, na.rm=T), y = max(tmp_wt$obs, na.rm=T),
           label = paste0("PCC: ", round(cor(tmp$obs, tmp$pred),2)),
           vjust=1, hjust=0, size=5)


pdf("Fig7A.pdf", height = 5.5, width = 10.5)
plot_grid(plotlist = list(gg_boxplot,gg_scater), nrow = 1, rel_widths = c(1,1.2), align = "h")
dev.off()



########
# Fig 7B
########

# function to plot logos of enhancers
library(ggseqlogo)
my.logo <- function(x, cutoff=NULL){
  p <- ggseqlogo(x, method='custom', seq_type='dna', ncol=1) +
    scale_x_continuous(breaks=seq(0,249,25), expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(panel.border = element_rect(colour="black", fill=NA),
          axis.ticks = element_line(colour="black"))
  if(!is.null(cutoff)) p <- p + geom_hline(yintercept = cutoff, lty="dashed")
  p
}

# load scores
twist_contr_scores <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/DeepSTARR_contr_scores_oligo_library.rds"))

pdf("Fig7B.pdf", height = 2.1, width = 13)
for(i in c("Synth_enh_dev_45", "Synth_enh_dev_53", "Synth_enh_dev_89")){
  p_dev <- my.logo(twist_contr_scores$dev[[grep(paste0(i,"_"), names(twist_contr_scores$dev))]]) + ggtitle(i) + scale_y_continuous("Dev scores")
  print(p_dev)
}
dev.off()
