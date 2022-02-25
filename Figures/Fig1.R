
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(GenomicRanges)
library(ggplot2)
theme_set(theme_light() + theme(axis.text = element_text(colour = "black")))
library(ggpointdensity)

########
# Fig 1D
########

df <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/STARRseq_data_DeepSTARR_predictions.txt")

# Prediction performance per set
scater_list_test <- lapply(c(dev="dev",
                             hk="hk"), function(class){
                               scater_list2 <- lapply(unique(df$set), function(set){
                                 tmp <- df[df$set %in% set,]
                                 
                                 if(class=="dev"){
                                   tmp$obs <- tmp$Dev_log2_enrichment
                                   tmp$pred <- tmp$Predictions_dev
                                 }else{
                                   tmp$obs <- tmp$Hk_log2_enrichment
                                   tmp$pred <- tmp$Predictions_hk
                                 }
                                 
                                 gg <- ggplot(tmp, aes(obs, pred)) +
                                   geom_pointdensity(size=0.3) +
                                   scale_color_gradient(low = "grey70", high = "grey20") +
                                   guides(col=F) +
                                   scale_x_continuous(paste0(class, " enhancer activity [log2]"),
                                                      breaks = seq(-12,20,2)) +
                                   scale_y_continuous(paste0("Predicted ", class, " activity [log2]"),
                                                      breaks = seq(-12,20,2)) +
                                   geom_abline(slope = 1, intercept = 0, linetype="dashed", col="grey30") +
                                   theme_bw() +
                                   theme(panel.background = element_rect(fill="white",colour="white"), panel.grid = element_blank(), axis.line=element_line(colour="black"),
                                         axis.text=element_text(size=14, colour="black"),
                                         axis.title=element_text(size=16, colour="black"),
                                         plot.title = element_text(size=18, hjust = 0.5, colour="black"), plot.subtitle = element_text(size=14, hjust = 0.5)) +
                                   annotate("text",  x=min(tmp$obs, na.rm=T), y = max(tmp$pred, na.rm=T),
                                            label = paste0("PCC: ", round(cor(tmp$obs, tmp$pred),2)),
                                            vjust=1, hjust=0, size=5)
                                 
                                 return(gg)
                               })
                               return(scater_list2)
                             })

# plot test set
ggsave("Fig1D_dev.png", scater_list_test$dev[[3]], height = 4, width = 4.2, type="cairo")
ggsave("Fig1D_hk.png", scater_list_test$hk[[3]], height = 4, width = 4.2, type="cairo")


########
# Fig 1E
########

## test set with actual peak summits
Peaks <- list(dev=read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/DSCP_200bp_gw.UMI_cut_merged.peaks.txt", header = F),
              hk=read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/RpS12_200bp_gw.UMI_cut_merged.peaks.txt", header = F))

Peaks <- lapply(Peaks, function(x){
  gr <- makeGRangesFromDataFrame(x, seqnames.field = "V1", start.field = "V2", end.field = "V2", seqinfo = Dmelanogaster@seqinfo, keep.extra.columns = T)
  mcols(gr) <- mcols(gr)[,5:7]
  names(mcols(gr)) <- c("Enrch.", "Corr_enrch", "p_value")
  gr <- gr[gr$Corr_enrch>3]
  return(gr)
})

# resize to 249bp (move 1nt because of bed file 0-based)
Peaks <- lapply(Peaks, function(x) resize(shift(x,-1), fix="center", 249))

# only test set
set <- "Test"
tmp <- df[df$set %in% "Test" & df$class %in% c("positive_peaks"),]

# get enhancers
tmp <- tmp[paste0(tmp$seqnames, tmp$start) %in% c(paste0(Peaks$dev@seqnames, Peaks$dev@ranges@start),
                                                  paste0(Peaks$hk@seqnames, Peaks$hk@ranges@start)),]

# calculate log2 FC
tmp$obs_log2FC <- tmp$Dev_log2_enrichment-tmp$Hk_log2_enrichment
tmp$pred_log2FC <- tmp$Predictions_dev-tmp$Predictions_hk

scater_list2_test <- ggplot(tmp, aes(obs_log2FC, pred_log2FC, col=obs_log2FC)) +
  geom_point(size=0.6) +
  scale_color_gradient2("log2FC [obs]", low = "#0C89CA", mid = "grey70", midpoint = 0, high = "#EF4A24") +
  scale_x_continuous("log2FC dev vs hk [observed]",
                     # limits = c(min(c(tmp$obs,0)),max(tmp$obs)),
                     breaks = seq(-12,20,2)) +
  scale_y_continuous("log2FC dev vs hk [predicted]",
                     # limits = c(min(c(tmp$pred,0)),max(tmp$pred)),
                     breaks = seq(-12,20,2)) +
  geom_abline(slope = 1, intercept = 0, linetype="dashed", col="grey30") +
  theme_bw() +
  theme(panel.background = element_rect(fill="white",colour="white"), panel.grid = element_blank(), axis.line=element_line(colour="black"),
        axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=16, colour="black"),
        plot.title = element_text(size=18, hjust = 0.5, colour="black"), plot.subtitle = element_text(size=14, hjust = 0.5))

ggsave("Fig1E.png", scater_list2_test, height = 4, width = 5.4, type="cairo")
