
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(ggplot2)
theme_set(theme_light() + theme(axis.text = element_text(colour = "black")))
library(patchwork)
library(dplyr)

motif_colours <- c("4 controls"="grey60",
                   "neg_ctrl_AAGTTT"="grey60",
                   "neg_ctrl_ATTAGT"="grey60",
                   "neg_ctrl_ATAAG"="grey60",
                   "neg_ctrl_AGTTG"="grey60",
                   "ctrl_AAGTTT"="grey60",
                   "ctrl_ATTAGT"="grey60",
                   "ctrl_ATAAG"="grey60",
                   "ctrl_AGTTG"="grey60",
                   "AP-1"=rgb(177,53,115, maxColorValue = 255),
                   "CREB1"="#377EB8",
                   "MAF"="#FF7F00", # 984EA3
                   "MECP2"="#984EA3", # E41A1C
                   "P53"="#4DAF4A",
                   "Ebox/MYC"="#FFFF33",
                   "EGR1"="#A65628",
                   "ETS"=rgb(255,102,102, maxColorValue = 255), # "#FF7F00"
                   "E2F1"="#F781BF",
                   "Random"="grey60")

########
# Fig 6B
########

p="5e-04"

# read table treated with PWM confident instances
Motif_instances_df <- readRDS(uel("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Hhuman_Motif_instances_df_treated_PWM_info_2.rds"))
Motif_instances_df <- Motif_instances_df[Motif_instances_df$Motif_2 %in% c("Random", "AP1", "P53", "MAF", "CREB1", "ETS", "EGR1"),]
Motif_instances_df$Motif_2 <- factor(Motif_instances_df$Motif_2,
                                     levels=c("Random", "AP1", "P53", "MAF", "CREB1", "ETS", "EGR1"),
                                     labels=c("4 controls", "AP-1", "P53", "MAF", "CREB1", "ETS", "EGR1"))

# Only enhancers that require the motifs (log2FC <= -1) and have a confident motif instance
# and neg control motifs
Motif_instances_df <- Motif_instances_df[(complete.cases(Motif_instances_df$log2FC_all_instances) & Motif_instances_df$log2FC_all_instances <= -1 & Motif_instances_df[,paste0("Motif_p", p)]>0) | Motif_instances_df$Motif_2 %in% "4 controls",]
table(Motif_instances_df$Motif_2, Motif_instances_df$`Motif_p5e-04`>0)
table(Motif_instances_df$Motif_2, Motif_instances_df$`Motif_p1e-04`>0)
table(Motif_instances_df$Motif_2, Motif_instances_df$Number_of_instances)

# calculate delta log2FC
df_delta <- Motif_instances_df %>%
  group_by(Sequence_ID, Motif) %>%
  summarise(Motif_2=unique(Motif_2),
            Number_of_instances=unique(Number_of_instances),
            Delta=max(log2FC)-min(log2FC),
            motif_counts=n())

table(df_delta$Motif, df_delta$Number_of_instances)
table(df_delta$Motif, df_delta$motif_counts)

# Plots

my_theme <- theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", size=15, angle=45, hjust=1),
        axis.text.y = element_text(colour="black", size=14),
        axis.title = element_text(colour="black", size=16),
        strip.text = element_text(colour="black", size=14))

# compare delta in FC to wt
tmp <- df_delta[df_delta$motif_counts>1,] # enhancers with at least 2 confident instances with measured mutation log2FC
tmp$Motif_2 <- droplevels(tmp$Motif_2)

# x labels counts per boxplot
xlabels <- sapply(as.character(levels(tmp$Motif_2)), function(x){
  paste0(gsub("_", " ", x)," (", table(tmp$Motif_2)[x], ")")
})

gg <- ggplot(tmp, aes(Motif_2, Delta, fill=Motif_2)) +
  geom_violin(alpha=0.7) +
  geom_boxplot(outlier.size = -1, color="black", width=0.15, size=0.8, fill=NA) +
  scale_fill_manual(values=motif_colours) +
  guides(fill=F) +
  scale_x_discrete("Mutated motifs (# Enhancers)", labels=xlabels) +
  scale_y_continuous(paste0("log2 FC between instances\nin the same enhancer"), breaks = seq(-10,10,1), expand=c(0,0), limits = c(0,5.5)) +
  my_theme

pdf("Fig6B.pdf", height = 5, width = 6)
print(gg + geom_hline(yintercept = 1, col="black", linetype="dashed", size=0.7))
dev.off()


########
# Fig 6C
########

tmp_melt <- reshape::melt(table(tmp$Motif_2, tmp$Delta >= log2(2)))
tmp_melt$Var.1 <- factor(tmp_melt$Var.1,
                         levels = c("4 controls", "AP-1", "P53", "MAF", "CREB1", "ETS", "EGR1"))

# x labels counts per boxplot
xlabels <- sapply(as.character(levels(tmp$Motif_2)), function(x){
  paste0(gsub("_", " ", x)," (", table(tmp$Motif_2)[x], ")")
})

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

av <- round(mean(prop.table(table(tmp$Motif_2, tmp$Delta >= log2(2)), margin = 1)[-1,2]),2)
ggsave("Fig6C.pdf",
       gg_bar + geom_bar(stat = "identity", position="fill", aes(fill = factor(Var.2)), width=0.8, colour="black", size=0.3) + geom_hline(yintercept = av, linetype="dashed", col="black") + scale_y_continuous("Proportion of enhancers [%]", expand=c(0,0), breaks = c(seq(0,1,0.25),av), labels = c(seq(0,1,0.25),av)*100),
       width = 6, height = 5)



########
# Fig 6D
########


### data from https://www.vierstra.org/resources/dgf
Motif_instances_gr <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Human_Motif_instances_gr_DNase_footprinting.rds"))

# compare log2FC and footprint in accessible enhancers
Motif_instances_gr_access <- Motif_instances_gr[Motif_instances_gr$DHS %in% "TRUE" & Motif_instances_gr$name %in% c("AP1", "P53", "MAF", "CREB1", "ETS", "EGR1")]

Motif_instances_gr_access$name <- factor(Motif_instances_gr_access$name,
                                         levels = c("ctrl_AAGTTT", "ctrl_ATTAGT", "ctrl_ATAAG", "ctrl_AGTTG",
                                                    "AP1", "P53", "MAF", "CREB1", "ETS", "EGR1"))

Motif_instances_gr_access$Motif_mutated2 <- as.character(Motif_instances_gr_access$name)
Motif_instances_gr_access$Motif_mutated2[grep("ctrl", Motif_instances_gr_access$Motif_mutated2)] <- "4 controls"
Motif_instances_gr_access$Motif_mutated2 <- factor(Motif_instances_gr_access$Motif_mutated2,
                                                   levels = c("AP1", "P53", "MAF", "CREB1", "ETS", "EGR1"),
                                                   labels = c("AP-1", "P53", "MAF", "CREB1", "ETS", "EGR1"))

# cell line
cell_line="h_RKO_DS40362"
# p-value of footprints
p="0.001"

tmp <- data.frame(Motif_instances_gr_access)
tmp$class <- tmp[,paste0(cell_line, "_fpr", p)]

# significance
pv <- sapply(levels(tmp$Motif_mutated2), function(TF) wilcox.test(tmp$score[tmp$class == "FALSE" & tmp$Motif_mutated2 %in% TF], tmp$score[tmp$class == "TRUE" & tmp$Motif_mutated2 %in% TF])$p.value)
names(pv) <- levels(tmp$Motif_mutated2)
print(paste0(cell_line, "_", p))
print(p.adjust(pv, "fdr"))

gg <- ggplot(tmp, aes(class, score, fill=Motif_mutated2, alpha=class)) +
  geom_boxplot(outlier.size = -1) +
  facet_grid(~Motif_mutated2) +
  scale_fill_manual(values=motif_colours) +
  scale_alpha_manual(values=c(0.5,1)) +
  guides(fill=F, alpha=F) +
  geom_hline(yintercept = 0, linetype="dashed", col="grey40") +
  scale_x_discrete("TF footprint", labels=c("-", "+")) +
  scale_y_continuous("log2FC enhancer activity", breaks = seq(-10,4,1)) +
  ggtitle(paste0("TF footprint in ", sapply(strsplit(cell_line,"_"), `[`, 2), " cells (fpr ", p, ")")) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(colour="black", size=15),
        axis.text.y = element_text(colour="black", size=14),
        axis.title = element_text(colour="black", size=16),
        strip.text = element_text(colour="black", size=13))

pdf("Fig6D.pdf", width = 7, height = 4)
print(gg)
dev.off()


########
# Fig 6E
########

### create linear models
library(caret)

# load main table
mutation_data_3 <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Human_mutation_data_corrected_instance_strand.rds"))
mutation_data_3$instance_center <- mutation_data_3$instance_start+(mutation_data_3$instance_end-mutation_data_3$instance_start)/2

#### extend motif core sequence in flanks to have enough nucleotides for motif matching
enh_gr <- GRanges(paste0(sapply(strsplit(mutation_data_3$Sequence_ID,"_"), `[`, 1),":",
                         sapply(strsplit(mutation_data_3$Sequence_ID,"_"), `[`, 2),"-",
                         sapply(strsplit(mutation_data_3$Sequence_ID,"_"), `[`, 3)),
                  strand = mutation_data_3$Strand)
motif_flank_length <- 5 # add 5bp on each side
mutation_data_3$instance_sequence_extended_5bp <- substr(as.character(getSeq(Hsapiens,enh_gr)),
                                                         mutation_data_3$instance_start-(motif_flank_length),
                                                         mutation_data_3$instance_end+(motif_flank_length))
mutation_data_3$instance_sequence_extended_5bp[mutation_data_3$instance_strand=="-"] <- as.character(reverseComplement(DNAStringSet(mutation_data_3$instance_sequence_extended_5bp[mutation_data_3$instance_strand=="-"])))


# motif positions
motif_scores_list_positions <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Human_Motif_oligo_positions_core_list.rds"))
motif_scores_list_positions <- do.call(rbind, motif_scores_list_positions)
motif_scores_list_positions$instance_center <- motif_scores_list_positions$start+(motif_scores_list_positions$end-motif_scores_list_positions$start)/2

# vars
vars <- c("Residuals",
          "Number_of_instances",
          "core",
          "flank_",
          "motif_position",
          "Close_FOSL1",
          "Close_TP53",
          "Close_MAFK",
          "Close_CREB1",
          "Close_ETS1",
          "Close_EGR1")
names(vars) <- c("Residuals",
                 "Number of instances",
                 paste0("Motif core"),
                 "Motif flanks",
                 "Motif position",
                 "Dist to AP-1",
                 "Dist to P53",
                 "Dist to MAF",
                 "Dist to CREB1",
                 "Dist to ETS",
                 "Dist to EGR1"
)

col_vars <- c("Residuals"="grey90",
              "Number of instances"="grey70",
              "Motif core"="grey50",
              "Motif flanks"="grey35",
              "Motif position"="#999900",
              "Dist to AP-1"=as.character(motif_colours["AP-1"]),
              "Dist to P53"=as.character(motif_colours["P53"]),
              "Dist to MAF"=as.character(motif_colours["MAF"]),
              "Dist to CREB1"=as.character(motif_colours["CREB1"]),
              "Dist to ETS"=as.character(motif_colours["ETS"]),
              "Dist to EGR1"=as.character(motif_colours["EGR1"])
)

m1_list <- as.character(unique(mutation_data_3$Motif_mutated)[-c(9:12)])
names(m1_list) <- m1_list

List_motif_model_results <- lapply(m1_list[c(1,3:6,9)], function(motif){
  
  ### get motif affinity
  if(motif=="CREB1") m <- "CREB1_MA0018.3"
  if(motif=="E2F1") m <- "E2F1_HUMAN.H11MO.0.A"
  if(motif=="EGR1") m <- "EGR1_C2H2_1"
  if(motif=="ETS") m <- "ETS1_HUMAN.H11MO.0.A"
  if(motif=="AP-1") m <- "FOSL1_MA0477.1"
  if(motif=="MAF") m <- "MAFK_MA0496.2"
  if(motif=="MECP2") m <- "MECP2_HUMAN.H11MO.0.C"
  if(motif=="Ebox/MYC") m <- "MYC_MA0147.3"
  if(motif=="P53") m <- "TP53_MA0106.3"
  
  df <- mutation_data_3[mutation_data_3$Motif_mutated %in% motif,c(1:3,12, 25, 28:33,35,36,44,45)]
  df$PWM <- mutation_data_3[mutation_data_3$Motif_mutated %in% motif,m]
  # df <- df[complete.cases(df$PWM),] # only instances with enough (10?) nucleotides in flanks
  
  # only instances with full flanks
  df <- df[nchar(df$instance_sequence_extended_5bp)==median(nchar(df$instance_sequence_extended_5bp)),]
  
  # split instances per column
  df2 <- as.data.frame(sapply(1:max(nchar(df$instance_sequence_extended_5bp)), function(i) sapply(strsplit(df$instance_sequence_extended_5bp,""), `[`, i)))
  names(df2) <- 1:ncol(df2)
  names(df2)[1:5] <- paste0("flank_left_", 5:1)
  names(df2)[names(df2) %in% 6:names(df2)[(ncol(df2)-4)]] <- paste0("core_", 1:length(6:names(df2)[(ncol(df2)-4)]))
  names(df2)[(ncol(df2)-4):ncol(df2)] <- paste0("flank_right_", 1:5)
  
  # remove positions with only 1 level
  df2 <- df2[,!sapply(1:ncol(df2), function(i) length(levels(df2[,i]))) == 1]
  
  df_final <- cbind(df, df2)
  
  ### get distance to closest instance of each motif type and bin distances
  for(m2 in c("FOSL1_MA0477.1", "TP53_MA0106.3", "CREB1_MA0018.3", "MAFK_MA0496.2", "EGR1_C2H2_1", "ETS1_HUMAN.H11MO.0.A")){
    df_final$dist_tmp <- sapply(1:nrow(df_final), function(i){
      dist <- df_final$instance_center[i]-motif_scores_list_positions$instance_center[motif_scores_list_positions$Sequence %in% df_final$Sequence_ID[i] & motif_scores_list_positions$Motif %in% m2]
      dist <- dist[!dist==0]
      out <- dist[abs(dist)==min(abs(dist))][1]
      if(length(out)>0){
        return(abs(out))
      }else{
        return(NA)}
    })
    df_final[,paste0("Close_", m2)] <- "NO"
    df_final[,paste0("Close_", m2)][which(df_final$dist_tmp < 25)] <- "Close"
    df_final[,paste0("Close_", m2)][which(df_final$dist_tmp >= 25 & df_final$dist_tmp <= 50)] <- "Intermediate"
    df_final[,paste0("Close_", m2)][which(df_final$dist_tmp > 50)] <- "Distal"
    df_final[,paste0("Close_", m2)] <- relevel(factor(df_final[,paste0("Close_", m2)]),
                                               "NO")
    # it is fine to keep dist_tmp column in table
  }
  
  # add position of the motif
  df_final$motif_position <- "center"
  df_final$motif_position[df_final$instance_end<=100 | df_final$instance_start>=150] <- "flanks" # outside middle 100bp
  df_final$motif_position[df_final$instance_end<=74 | df_final$instance_start>=175] <- "boundaries" # outside middle 100bp
  df_final$motif_position <- factor(df_final$motif_position, levels = c("flanks", "center", "boundaries"))
  table(df_final$motif_position)
  
  df_final$instance_center <- df_final$instance_start+(df_final$instance_end-df_final$instance_start)/2
  df_final$motif_position_bins <- cut(125-df_final$instance_center,breaks = seq(-125,125,50))
  
  # multiple models (with caret) with motif counts, presence of other motifs, and motif flanks
  # 10-fold CV
  trctrl <- trainControl(method = "repeatedcv",
                         number = 10,
                         repeats = 1,
                         classProbs = T,
                         savePredictions = T)
  model_list <- c("lm")
  names(model_list) <- model_list
  ### model with motif distance and NO
  model_predictions_dist_NO <- lapply(model_list, function(model){
    set.seed(1234)
    model_trained <- train(log2FC_wt_mut~., data=df_final[complete.cases(df_final$log2FC_wt_mut),grep("Number_of|Close_|flank|core|motif_position$|^log2FC_wt_mut", names(df_final))], method = model,
                           trControl=trctrl,
                           importance = TRUE,
                           tuneLength = 1)
    pred <- model_trained$pred
    return(list(Model_trained=model_trained,
                Pred=pred,
                adj.r.squared=summary(lm(pred$pred~pred$obs))$adj.r.squared,
                PCC=cor.test(pred$pred, pred$obs)$estimate))
  })
  
  ### 100 models - to get feature importance
  df1 <- df_final[complete.cases(df_final$log2FC_wt_mut),grep("Number_of|Close_|flank|core|motif_position|^log2FC_wt_mut", names(df_final))]
  var_expl_all <- sapply(1:100, function(i){
    df3 <- df1[,sample(ncol(df1))]
    fit <- lm(log2FC_wt_mut~., data=df3)
    af <- anova(fit)
    af$PctExp <- af$"Sum Sq"/sum(af$"Sum Sq")*100
    var_expl <- sapply(vars, function(x) sum(af$PctExp[grep(x, rownames(af))]))
    var_expl
  })
  var_expl_all <- rowMeans(var_expl_all)
  
  
  ### save
  return(list(Data=df_final,
              Models_dist_NO=model_predictions_dist_NO,
              adj.r.squared_dist_NO=sapply(model_predictions_dist_NO, function(x) x$adj.r.squared),
              PCC_dist_NO=sapply(model_predictions_dist_NO, function(x) x$PCC),
              lm_var_expl_average=var_expl_all,
              boxplot_position=gg_position
  ))
  
})
saveRDS(List_motif_model_results, "List_motif_model_results_flank_core_motif_position_sequences.rds")


### summarise importance of different features per motif type - for linear model - easier to interpret
List_motif_model_results <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Human_List_motif_model_results_flank_core_motif_position_sequences.rds"))

vars <- c("Number_of_instances",
          "flank_left_",
          "core",
          "flank_right_",
          "motif_position",
          "Close_FOSL1",
          "Close_TP53",
          "Close_MAFK",
          "Close_CREB1",
          "Close_ETS1",
          "Close_EGR1")
names(vars) <- c("Number of instances",
                 "Flanking left",
                 paste0("Motif core"),
                 "Flanking right",
                 "Motif position",
                 "Dist to AP-1",
                 "Dist to P53",
                 "Dist to MAF",
                 "Dist to CREB1",
                 "Dist to ETS",
                 "Dist to EGR1"
)

motif_anno <- data.frame(row.names = names(sapply(List_motif_model_results, function(m) m$adj.r.squared_dist_NO["lm"])),
                         # Model.adj.R2 = sapply(List_motif_model_results, function(m) m$adj.r.squared_dist_NO["lm"]),
                         PCC = sapply(List_motif_model_results, function(m) m$PCC_dist_NO["lm.cor"]))
rownames(motif_anno) <- gsub(".lm", "", rownames(motif_anno))

# Get FDR-corrected p-value bin
summary_lm <- do.call(rbind, lapply(names(List_motif_model_results), function(model){
  # Linear model - significance
  out <- as.data.frame(summary(List_motif_model_results[[model]]$Models_dist_NO$lm$Model_trained)$coefficients)
  out[,"Motif"] <- model
  out
}))
# summary_lm$FDR <- p.adjust(summary_lm$`Pr(>|t|)`, method="fdr")

summary_lm <- sapply(unique(summary_lm$Motif), function(model){
  
  # Linear model - significance
  m <-summary_lm[summary_lm$Motif %in% model,]
  var_imp <- sapply(vars, function(x){
    p <- min(m[grep(x, rownames(m)),4])
    # if(p>=0.05) p2 <- 1
    # if(p<0.05) p2 <- 0.05
    # if(p<0.01) p2 <- 0.01
    # if(p<0.001) p2 <- 0.001
    # if(p<0.0001) p2 <- 0.0001
    if(p>=0.05) p2 <- ">0.05"
    if(p<0.05) p2 <- "0.05"
    if(p<0.01) p2 <- "0.01"
    if(p<0.001) p2 <- "0.001"
    if(p<0.0001) p2 <- "0.0001"
    p2
  })
  
  return(var_imp)
  
})
summary_lm


### heatmap
# motifs in rows, features in columns
# add PCC
out <- reshape2::melt(summary_lm)
out$Var2 <- factor(out$Var2, levels=rev(c("AP-1", "P53", "MAF", "CREB1", "ETS", "EGR1")))
out$value <- factor(out$value, levels=c(">0.05",
                                        "0.05",
                                        "0.01",
                                        "0.001",
                                        "0.0001"))
g <- ggplot(out, aes(Var1, Var2, fill=value)) +
  geom_tile(col="black") + 
  scale_fill_manual("linear model p-value", values=c(">0.05"="white",
                                                     "0.05"="#FDDBC7",
                                                     "0.01"="#F4A582",
                                                     "0.001"="#D6604D",
                                                     "0.0001"="#B2182B")) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "bottom")

motif_anno$Var2 <- factor(rownames(motif_anno), levels=rev(c("AP-1", "P53", "MAF", "CREB1", "ETS", "EGR1")))

g_PCC <- ggplot(motif_anno, aes("1", Var2, fill=PCC)) +
  geom_tile(col="black") + 
  scale_fill_gradient("PCC", low = "white", high = "green4", breaks = seq(0,10,0.2), limits=c(0,0.7)) +
  guides(fill=F) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank())

ggsave("Fig6E.pdf",
       g_PCC + g + plot_layout(widths=c(1,10)),
       width = 7.5, height = 5)


#### plot for all features - raw importance
vars <- c("Number_of_instances",
          paste0("flank_left_5", c("C", "G", "T")),
          paste0("flank_left_4", c("C", "G", "T")),
          paste0("flank_left_3", c("C", "G", "T")),
          paste0("flank_left_2", c("C", "G", "T")),
          paste0("flank_left_1", c("C", "G", "T")),
          "core",
          paste0("flank_right_1", c("C", "G", "T")),
          paste0("flank_right_2", c("C", "G", "T")),
          paste0("flank_right_3", c("C", "G", "T")),
          paste0("flank_right_4", c("C", "G", "T")),
          paste0("flank_right_5", c("C", "G", "T")),
          paste0("motif_position", c("center", "boundaries")),
          paste0("Close_FOSL1_MA0477.1", c("Close", "Intermediate", "Distal")),
          paste0("Close_TP53_MA0106.3", c("Close", "Intermediate", "Distal")),
          paste0("Close_MAFK_MA0496.2", c("Close", "Intermediate", "Distal")),
          paste0("Close_CREB1_MA0018.3", c("Close", "Intermediate", "Distal")),
          paste0("Close_ETS1_HUMAN.H11MO.0.A", c("Close", "Intermediate", "Distal")),
          paste0("Close_EGR1_C2H2_1", c("Close", "Intermediate", "Distal")))
names(vars) <- c("Number of instances",
                 paste0("Flanking nt -5 ", c("C", "G", "T")),
                 paste0("Flanking nt -4 ", c("C", "G", "T")),
                 paste0("Flanking nt -3 ", c("C", "G", "T")),
                 paste0("Flanking nt -2 ", c("C", "G", "T")),
                 paste0("Flanking nt -1 ", c("C", "G", "T")),
                 paste0("Motif core"),
                 paste0("Flanking nt +1 ", c("C", "G", "T")),
                 paste0("Flanking nt +2 ", c("C", "G", "T")),
                 paste0("Flanking nt +3 ", c("C", "G", "T")),
                 paste0("Flanking nt +4 ", c("C", "G", "T")),
                 paste0("Flanking nt +5 ", c("C", "G", "T")),
                 paste0("Motif position - ", c("center", "boundaries")),
                 paste0("Dist to AP-1 - ", c("close", "intermediate", "distal")),
                 paste0("Dist to P53 - ", c("close", "intermediate", "distal")),
                 paste0("Dist to MAF - ", c("close", "intermediate", "distal")),
                 paste0("Dist to CREB1 - ", c("close", "intermediate", "distal")),
                 paste0("Dist to ETS - ", c("close", "intermediate", "distal")),
                 paste0("Dist to E2F1 - ", c("close", "intermediate", "distal"))
)

# Get FDR-corrected p-value bin
summary_lm <- do.call(rbind, lapply(names(List_motif_model_results), function(model){
  # Linear model - significance
  out <- as.data.frame(summary(List_motif_model_results[[model]]$Models_dist_NO$lm$Model_trained)$coefficients)
  out[,"Motif"] <- model
  out
}))
# summary_lm$FDR <- p.adjust(summary_lm$`Pr(>|t|)`, method="fdr")

summary_lm <- sapply(unique(summary_lm$Motif), function(model){
  
  # Linear model - significance
  m <-summary_lm[summary_lm$Motif %in% model,]
  var_imp <- sapply(vars, function(x){
    m2 <- m[grep(x, rownames(m)),]
    m2 <- m2[order(m2[,4]),]
    if(x=="core" & length(which(m2[m2$Estimate<0,6] < 0.1))>0) m2 <- m2[m2$Estimate<0,] # for core show positive associations if possible
    
    p <- m2[1,4]
    # if(p>=0.05) p2 <- 1
    # if(p<0.05) p2 <- 0.05
    # if(p<0.01) p2 <- 0.01
    # if(p<0.001) p2 <- 0.001
    # if(p<0.0001) p2 <- 0.0001
    if(p>=0.05) p2 <- ">0.05"
    if(p<0.05) p2 <- "0.05"
    if(p<0.01) p2 <- "0.01"
    if(p<0.001) p2 <- "0.001"
    if(p<0.0001) p2 <- "0.0001"
    
    if(sign(m2[1,1])==1 & p<0.1) p2 <- paste0(p2, " (-)") # revert signal, because I want positive contribution to correlate with stronger negative FC (of mutation)
    if(sign(m2[1,1])==-1 & p<0.1) p2 <- paste0(p2, " (+)")
    p2
  })
  
  return(var_imp)
  
})
summary_lm



### heatmap
# motifs in rows, features in columns
# add PCC
out <- reshape2::melt(summary_lm)
out$Var2 <- factor(out$Var2, levels=rev(c("AP-1", "P53", "MAF", "CREB1", "ETS", "EGR1")))
out$value <- factor(out$value, levels=c("0.0001 (-)",
                                        "0.001 (-)",
                                        "0.01 (-)",
                                        "0.05 (-)",
                                        ">0.05",
                                        "0.05 (+)",
                                        "0.01 (+)",
                                        "0.001 (+)",
                                        "0.0001 (+)"))

g <- ggplot(out, aes(Var1, Var2, fill=value)) +
  geom_tile(col="black") + 
  scale_fill_manual("linear model\np-value (direction)", values=c("0.0001 (-)"="#2166AC",
                                                                  "0.001 (-)"="#4393C3",
                                                                  "0.01 (-)"="#92C5DE",
                                                                  "0.05 (-)"="#D1E5F0",
                                                                  ">0.05"="white",
                                                                  "0.05 (+)"="#FDDBC7",
                                                                  "0.01 (+)"="#F4A582",
                                                                  "0.001 (+)"="#D6604D",
                                                                  "0.0001 (+)"="#B2182B"),
                    breaks=c("0.0001 (-)",
                             "0.001 (-)",
                             "0.01 (-)",
                             "0.05 (-)",
                             ">0.05",
                             "0.0001 (+)",
                             "0.001 (+)",
                             "0.01 (+)",
                             "0.05 (+)"
                    ), drop=F) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=10, angle=90, hjust=1),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=15),
        legend.position = "bottom")

motif_anno$Var2 <- factor(rownames(motif_anno), levels=rev(c("AP-1", "P53", "MAF", "CREB1", "ETS", "EGR1")))

g_PCC <- ggplot(motif_anno, aes("1", Var2, fill=PCC)) +
  geom_tile(col="black") + 
  scale_fill_gradient("PCC", low = "white", high = "green4", breaks = seq(0,10,0.2), limits=c(0,0.7)) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(colour="black", size=14),
        legend.title = element_text(colour="black", size=15)
  )

ggsave("FigS21B.pdf",
       g + g_PCC + plot_layout(widths=c(23,1)),
       width = 13, height = 6)


########
# Fig 6F
########

# load main table
mutation_data_3 <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Human_mutation_data_3.rds"))
mutation_data_3$instance_center <- mutation_data_3$instance_start+(mutation_data_3$instance_end-mutation_data_3$instance_start)/2

# motif positions
motif_scores_list_positions <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Human_Motif_oligo_positions_core_list.rds"))
motif_scores_list_positions <- do.call(rbind, motif_scores_list_positions)
motif_scores_list_positions$instance_center <- motif_scores_list_positions$start+(motif_scores_list_positions$end-motif_scores_list_positions$start)/2

m1_list <- as.character(unique(mutation_data_3$Motif_mutated)[-c(9:12)])
names(m1_list) <- m1_list

pdf("Fig6F.pdf", width = 3, height = 4.9)
for(motif in m1_list[c(5,4)]){
  
  df <- mutation_data_3[mutation_data_3$Motif_mutated %in% motif, c(1:3,12, 25, 28:32, 46)]
  
  ### get motif affinity
  if(motif=="CREB1") m <- "CREB1_MA0018.3"
  if(motif=="E2F1") m <- "E2F1_HUMAN.H11MO.0.A"
  if(motif=="EGR1") m <- "EGR1_C2H2_1"
  if(motif=="ETS") m <- "ETS1_HUMAN.H11MO.0.A"
  if(motif=="AP1") m <- "FOSL1_MA0477.1"
  if(motif=="MAF") m <- "MAFK_MA0496.2"
  if(motif=="MECP2") m <- "MECP2_HUMAN.H11MO.0.C"
  if(motif=="Ebox/MYC") m <- "MYC_MA0147.3"
  if(motif=="P53") m <- "TP53_MA0106.3"
  
  ### get distance to closest instance of each motif type
  
  m2_list <- unique(motif_scores_list_positions$Motif)
  for(m2 in m2_list){
    df[,paste0(m2, "_dist")] <- sapply(1:nrow(df), function(i){
      dist <- df$instance_center[i]-motif_scores_list_positions$instance_center[motif_scores_list_positions$Sequence %in% df$Sequence_ID[i] & motif_scores_list_positions$Motif %in% m2]
      dist <- dist[!dist==0]
      out <- dist[abs(dist)==min(abs(dist))][1]
      return(out)
    })
  }
  
  # pair with each motif type
  for(m2 in m2_list[c(1,8)]){
    
    gg_df <- df[complete.cases(df[,paste0(m2, "_dist")]),] # restricting by complete distances means it has at least 2 when is homotypic pair, and at least one of each when heterotypic. Should I limit homotypic to 2 instances
    
    # remove overlapping instances (distance between centers > motif length)
    gg_df <- gg_df[abs(gg_df[,paste0(m2, "_dist")]) > median(gg_df$instance_width),]
    
    gg_df$class <- NA
    gg_df$class[abs(gg_df[,paste0(m2, "_dist")])<25] <- "<25"
    gg_df$class[abs(gg_df[,paste0(m2, "_dist")])>50] <- ">50"
    table(gg_df$class)
    
    m1 <- motif
    if(motif=="AP1") m1 <- "AP-1"
    
    if(m2=="FOSL1_MA0477.1") m2_label <- "AP-1"
    if(m2=="ETS1_HUMAN.H11MO.0.A") m2_label <- "ETS"
    
    g_final <- ggplot(gg_df[complete.cases(gg_df$class),], aes(class, log2FC_wt_mut, colour=class)) +
      geom_boxplot(fill=NA, size=1, outlier.size = -1) +
      scale_colour_manual(values=c("<25"="#2166AC", ">50"="#B2182B")) +
      guides(col=F) +
      scale_y_continuous(paste0("log2 FC enhancer activity [",m1,"]"), breaks = seq(-10,10,1)) +
      scale_x_discrete(paste0(m1,"/", m2_label," distance (bp)"),
                       labels=c(paste0("<=25\n(n=", length(which(gg_df$class=="<25")),")"),
                                paste0(">50\n(n=", length(which(gg_df$class==">50")),")"))) +
      geom_hline(yintercept = 0, linetype="dashed", col="grey60") +
      ggtitle(paste0("wilcox p: ", signif(wilcox.test(gg_df$log2FC~gg_df$class)$p.value,2)))
    
    if(motif=="AP1" & m2 =="FOSL1_MA0477.1") g_final <- g_final + coord_cartesian(ylim = c(0.7, -4.2))
    if(motif=="AP1" & m2 =="ETS1_HUMAN.H11MO.0.A") g_final <- g_final + coord_cartesian(ylim = c(0.7, -4.2))
    if(motif=="ETS" & m2 =="FOSL1_MA0477.1") g_final <- g_final + coord_cartesian(ylim = c(1.2, -2.6))
    
    print(g_final)
    
  }
  
}
dev.off()

