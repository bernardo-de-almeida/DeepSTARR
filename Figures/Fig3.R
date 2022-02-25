
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(GenomicRanges)
library(ggplot2)
theme_set(theme_light() + theme(axis.text = element_text(colour = "black")))
library(ggpointdensity)

motif_colours <- c(GATA=rgb(116,159,242, maxColorValue = 255),
                   GATAA=rgb(116,159,242, maxColorValue = 255),
                   TGA.TCA=rgb(177,53,115, maxColorValue = 255),
                   AP1=rgb(177,53,115, maxColorValue = 255),
                   "AP-1"=rgb(177,53,115, maxColorValue = 255),
                   GAGA=rgb(212,147,91, maxColorValue = 255),
                   Trl=rgb(212,147,91, maxColorValue = 255),
                   twist=rgb(121,170,109, maxColorValue = 255),
                   ATCGAT=rgb(0,128,255, maxColorValue = 255),
                   Dref=rgb(0,128,255, maxColorValue = 255),
                   Random="grey60",
                   "3 controls"="grey60",
                   ctrl_TAGG="grey60",
                   ctrl_GGGCT="grey60",
                   ctrl_CCTTA="grey60",
                   GGGCT="grey60"
)

########
# Fig 3A
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

# nucleotide contribution scores for oligos
twist_contr_scores <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/DeepSTARR_contr_scores_oligo_library.rds"))

# plots
pdf("Fig3A.pdf", height = 2.1, width = 13)

i="chr3L_13015084_13015332_-_wt_dCP"
my.logo(twist_contr_scores$dev[[grep(i, names(twist_contr_scores$dev), fixed = T)]]) + ggtitle(paste0(i, " - dev scores"))

dev.off()


########
# Fig 3B
########

mutation_data_and_DeepSTARR <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_mutation_data_3_and_DeepSTARR.rds"))
mutation_data_and_DeepSTARR <- mutation_data_and_DeepSTARR[,c(1:14,60:69,75,80,78,73,79,85,86)]
names(mutation_data_and_DeepSTARR)[30:31] <- c("DeepSTARR_dev", "DeepSTARR_hk")

Motifs <- list(GATAA=data.frame(Motif="GATAA", ID="flyfactorsurvey__srp_SANGER_5_FBgn0003507"),
               TGA.TCA=data.frame(Motif="TGA.TCA", ID="jaspar__MA0476.1"),
               twist=data.frame(Motif="twist", ID="flyfactorsurvey__twi_da_SANGER_5_FBgn0000413"),
               Trl=data.frame(Motif="Trl", ID="flyfactorsurvey__Trl_FlyReg_FBgn0013263"),
               ATCGAT=data.frame(Motif="ATCGAT", ID="homer__AVYTATCGATAD_DREF"))
Motifs <- do.call(rbind, Motifs)

pdf("Fig3B.pdf", width = 6.7, height = 4.7)
library(patchwork)
summary_statistics <- data.frame()
for(motif in unique(Motifs$Motif)){
  
  class <- ifelse(motif=="ATCGAT", "hk", "dev")
  high_col <- motif_colours[motif]
  
  if(motif=="twist"){
    tmp <- rbind(mutation_data_and_DeepSTARR[mutation_data_and_DeepSTARR$enhancer_group %in% class & mutation_data_and_DeepSTARR$Motif_mutated %in% "CA..TG" & mutation_data_and_DeepSTARR$wt_instance %in% c("CATCTG", "CAGATG", "CATATG"),],
                 mutation_data_and_DeepSTARR[mutation_data_and_DeepSTARR$enhancer_group %in% class & mutation_data_and_DeepSTARR$Motif_mutated %in% "twist",])
  }else if(motif=="Trl"){
    tmp <- mutation_data_and_DeepSTARR[mutation_data_and_DeepSTARR$enhancer_group %in% class & mutation_data_and_DeepSTARR$Motif_mutated %in% motif,]
    tmp <- tmp[c(intersect(grep("GAGAG", tmp$instance_sequence_extended), grep("+", tmp$instance_strand)),
                 intersect(grep("CTCTC", tmp$instance_sequence_extended), grep("-", tmp$instance_strand))),]
  }else{
    tmp <- mutation_data_and_DeepSTARR[mutation_data_and_DeepSTARR$enhancer_group %in% class & mutation_data_and_DeepSTARR$Motif_mutated %in% motif,]
  }
  
  if(motif=="ATCGAT"){model_list <- c(names(tmp)[31], as.character(Motifs$ID[Motifs$Motif %in% motif]))}else{model_list <- c(names(tmp)[30], as.character(Motifs$ID[Motifs$Motif %in% motif]))}
  
  for(model in model_list){
    
    if(class=="dev") tmp$var <- tmp$dev_log2FC_wt_mut
    if(class=="hk") tmp$var <- tmp$hk_log2FC_wt_mut
    tmp <- tmp[complete.cases(tmp$var),]
    
    g1 <- ggplot(tmp, aes("1", var)) +
      geom_violin(alpha=0.7, fill=high_col) +
      geom_boxplot(outlier.size = -1, color="black", width=0.18, size=0.8, fill=NA) +
      scale_y_continuous(paste0("log2 FC - ", motif, " mutant [observed]"), breaks = seq(-10,4,2)) +
      xlab(paste0("All instances\n(n=", nrow(tmp[complete.cases(tmp$var),]), ")")) +
      geom_hline(yintercept = 0, linetype="dashed", col="grey40") +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank())
    
    pc <- cor.test(tmp[,model], tmp$var,
                   method = "pearson", use="complete.obs")
    
    g2 <- ggplot(tmp, aes(tmp[,model], var)) +
      geom_pointdensity() +
      scale_color_gradient(low = "grey70", high = high_col) +
      guides(col=F) +
      scale_y_continuous(paste0("log2 FC - ", motif, " mutant"), breaks = seq(-10,4,2)) +
      scale_x_continuous(paste0("Predicted motif importance")) +
      ggtitle(paste0(motif, " mutations - ", model)) +
      geom_hline(yintercept = 0, linetype="dashed", col="grey40") +
      theme(axis.title.y = element_blank()) +
      annotate("text",  x=min(tmp[,model], na.rm=T), y = max(tmp$var, na.rm=T), label = paste0("PCC: ", round(pc$estimate,2)), vjust=1, hjust=0, size=5)
    
    print(g1 + g2 + plot_layout(widths = c(1,5)))
    
    summary_statistics <- rbind(summary_statistics, data.frame(motif=motif,
                                                               model=model,
                                                               PCC=pc$estimate,
                                                               pvalue=pc$p.value))
  }
}

dev.off()


########
# Fig 3C
########

# read table treated with PWM confident instances
Motif_instances_df <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_Motif_instances_df_treated_PWM_info_2.rds"))

table(Motif_instances_df$Motif, Motif_instances_df$`Motif_p5e-04`>0)
table(Motif_instances_df$Motif_2, Motif_instances_df$`Motif_p5e-04`>0)

# Only enhancers that require the motifs (log2FC <= -1) and have a confident motif instance
# Motif_instances_df <- Motif_instances_df[complete.cases(Motif_instances_df$log2FC_all_instances) & Motif_instances_df$log2FC_all_instances <= -1 & Motif_instances_df$`Motif_p5e-04`>0,]
# don't subset on PWM score and include  neg control motifs
Motif_instances_df <- Motif_instances_df[(complete.cases(Motif_instances_df$log2FC_all_instances) & Motif_instances_df$log2FC_all_instances <= -1) | Motif_instances_df$Motif_2 %in% "3 controls",]
table(Motif_instances_df$Motif_2, Motif_instances_df$`Motif_p5e-04`>0)

# calculate delta log2FC
df_delta <- Motif_instances_df %>%
  group_by(Sequence_ID, Motif) %>%
  summarise(Motif_2=unique(Motif_2),
            enhancer_group=unique(enhancer_group),
            Delta=max(log2FC)-min(log2FC),
            motif_counts=n())

table(df_delta$Motif, df_delta$motif_counts)
table(df_delta$Motif_2, df_delta$motif_counts)

my_theme <- theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", size=15, angle=45, hjust=1),
        axis.text.y = element_text(colour="black", size=14),
        axis.title = element_text(colour="black", size=16),
        strip.text = element_text(colour="black", size=14))

#### compare delta in FC to wt

tmp <- df_delta[df_delta$Motif_2 %in% c("3 controls", "AP-1", "GATA", "twist", "Trl", "Dref") & df_delta$motif_counts>1,]
tmp$Motif_2 <- droplevels(tmp$Motif_2)

# x labels counts per boxplot
xlabels <- sapply(as.character(levels(tmp$Motif_2)), function(x){
  paste0(gsub("_", " ", x)," (", table(tmp$Motif_2)[x], ")")
})

# test difference to 3 controls
wilcox.test(tmp$Delta[tmp$Motif_2 == "AP-1" & tmp$enhancer_group=="dev"], tmp$Delta[tmp$Motif_2 == "3 controls" & tmp$enhancer_group=="dev"]) # n.s.
wilcox.test(tmp$Delta[tmp$Motif_2 == "GATA" & tmp$enhancer_group=="dev"], tmp$Delta[tmp$Motif_2 == "3 controls" & tmp$enhancer_group=="dev"]) # ***
wilcox.test(tmp$Delta[tmp$Motif_2 == "twist" & tmp$enhancer_group=="dev"], tmp$Delta[tmp$Motif_2 == "3 controls" & tmp$enhancer_group=="dev"]) # *
wilcox.test(tmp$Delta[tmp$Motif_2 == "Trl" & tmp$enhancer_group=="dev"], tmp$Delta[tmp$Motif_2 == "3 controls" & tmp$enhancer_group=="dev"]) # ***
wilcox.test(tmp$Delta[tmp$Motif_2 == "Dref" & tmp$enhancer_group=="hk"], tmp$Delta[tmp$Motif_2 == "3 controls" & tmp$enhancer_group=="hk"]) # ****

gg <- ggplot(tmp, aes(Motif_2, Delta, fill=Motif_2)) +
  geom_violin(alpha=0.7) +
  geom_boxplot(outlier.size = -1, color="black", width=0.15, size=0.8, fill=NA) +
  facet_grid(~enhancer_group, scales = "free_x", space="free_x") +
  scale_fill_manual(values=motif_colours) +
  guides(fill=F) +
  scale_x_discrete("Mutated motifs", labels=xlabels) +
  scale_y_continuous(paste0("log2 FC between instances\nin the same enhancer"), breaks = seq(-10,10,1), expand=c(0,0), limits = c(0,5.9)) +
  my_theme

pdf(paste0("Fig3C_left.pdf"), height = 5, width = 7)
print(gg + geom_hline(yintercept = 1, col="black", linetype="dashed", size=0.7))
dev.off()

### barplots
tmp_melt <- reshape::melt(table(tmp$Motif_2, tmp$Delta >= log2(2)))
tmp_melt <- tmp_melt[tmp_melt$Var.1 %in% c("3 controls", "AP-1", "GATA", "twist", "Trl", "Dref"),]
tmp_melt$Var.1 <- factor(tmp_melt$Var.1, levels=c("3 controls", "AP-1", "GATA", "twist", "Trl", "Dref"))

gg_bar <- ggplot(tmp_melt, aes(Var.1, value)) +
  xlab("Mutated motifs") +
  scale_x_discrete("Mutated motifs", labels=xlabels) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour=motif_colours[levels(tmp_melt$Var.1)], size=15, angle=45, hjust=1),
        axis.text.y = element_text(colour="black", size=14),
        axis.title = element_text(colour="black", size=16),
        strip.text = element_text(colour="black", size=14)) +
  scale_fill_manual("> 2-fold diff", values=c("grey90", "grey50")) +
  scale_colour_manual("log2 FC < -1", values=c("grey90", "red")) +
  guides(colour=F)

# barplots - with line on average (excluding neg regions)
av <- round(mean(prop.table(table(tmp$Motif_2, tmp$Delta >= log2(2)), margin = 1)[-1,2]),3)
ggsave("Fig3C_right.pdf",
       gg_bar + geom_bar(stat = "identity", position="fill", aes(fill = factor(Var.2)), width=0.8, colour="black", size=0.3) + geom_hline(yintercept = av, linetype="dashed", col="black") + scale_y_continuous("Proportion of enhancers [%]", expand=c(0,0), breaks = c(seq(0,1,0.25),av), labels = c(seq(0,1,0.25),av)*100),
       width = 6, height = 5)


########
# Fig 3D
########

final_statistics <- summary_statistics
final_statistics$motif <- factor(final_statistics$motif, levels = c("GATAA", "TGA.TCA", "twist", "Trl", "ATCGAT"),
                                 labels = c("GATA", "AP-1", "twist", "Trl", "Dref"))

final_statistics$model2 <- "PWM"
final_statistics$model2[grep("DeepSTARR_dev", final_statistics$model)] <- "DeepSTARR log2FC (dev)"
final_statistics$model2[grep("DeepSTARR_hk", final_statistics$model)] <- "DeepSTARR log2FC (hk)"
final_statistics$model2 <- factor(final_statistics$model2, levels = c("PWM", "DeepSTARR log2FC (dev)", "DeepSTARR log2FC (hk)"))

# correct PCC from log2FC comparisons
final_statistics$PCC[-grep("DeepSTARR", final_statistics$model)] <- -final_statistics$PCC[-grep("DeepSTARR", final_statistics$model)]

gg <- ggplot(final_statistics, aes(x=motif, y=PCC, fill=model2)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.7, colour="black", size=0.2) +
  scale_fill_manual("Model", values=c("grey60", "orangered", "dodgerblue")) +
  scale_y_continuous("PCC", breaks = seq(-10,4,0.1), limits = c(0,0.55), expand = c(0,0)) +
  scale_x_discrete("Motif mutated", labels=c("GATA", "AP-1", "twist", "Trl", "Dref"))

ggsave("Fig3D.pdf", gg, width = 7.1, height = 4)
