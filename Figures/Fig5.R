
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(GenomicRanges)
library(ggplot2)
theme_set(theme_light() + theme(axis.text = element_text(colour = "black")))
library(patchwork)
library(dplyr)

high_col_list <- c(GATA=rgb(116,159,242, maxColorValue = 255),
                   GATAA=rgb(116,159,242, maxColorValue = 255),
                   TGA.TCA=rgb(177,53,115, maxColorValue = 255),
                   GAGA=rgb(212,147,91, maxColorValue = 255),
                   CATCTG=rgb(121,170,109, maxColorValue = 255),
                   CA..TG=rgb(121,170,109, maxColorValue = 255),
                   ATCGAT=rgb(0,128,255, maxColorValue = 255),
                   "AP-1"=rgb(177,53,115, maxColorValue = 255),
                   AP1=rgb(177,53,115, maxColorValue = 255),
                   Trl=rgb(212,147,91, maxColorValue = 255),
                   twist=rgb(121,170,109, maxColorValue = 255),
                   ETS=rgb(255,102,102, maxColorValue = 255),
                   SREBP=rgb(102,102,0, maxColorValue = 255),
                   Dref=rgb(0,128,255, maxColorValue = 255),
                   Ohler1=rgb(0,76,153, maxColorValue = 255),
                   Ohler5="#4C0099",
                   "Ebox/Ohler5"="#4C0099",
                   Ohler6=rgb(0,153,153, maxColorValue = 255),
                   Ohler7="#C899F7",
                   Random="grey60",
                   GGGCT="grey60")

########
# Fig 5B
########

### load motif pair co-occurence results
fisher_results <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_motif_pair_distance_preferences_genomic_enhancers_Fisher_enrichment.rds"))
fisher_results <- fisher_results %>% group_by(type) %>% mutate(fisher_FDR=p.adjust(fisher_p, method="fdr"),
                                                               fisher_FDR_label=ifelse(p.adjust(fisher_p, method="fdr")<0.05, "*", ""))
fisher_results <- fisher_results[fisher_results$type %in% "All_pairs",]
fisher_results$Distance <- factor(fisher_results$Distance, levels = c("0-25", "25-50", "50-75", "75-100", "100-125", "125-150", "150-250"))
fisher_results$m1 <- gsub("AP1", "AP-1", fisher_results$m1)
fisher_results$m2 <- gsub("AP1", "AP-1", fisher_results$m2)


### load linear regression association with enhancer activity
lm_dist <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_validation_activity_motif_distance_genomic_sequences_all_per_bin.txt", stringsAsFactors = F)
lm_dist <- lm_dist %>% mutate(FDR=p.adjust(pvalue, method="fdr"),
                              FDR_label=ifelse(p.adjust(pvalue, method="fdr")<0.05, "*", ""))

lm_dist$cutoff <- factor(lm_dist$cutoff, levels = c("0-25", "25-50", "50-75", "75-100", "100-125", "125-150", "150-250"))
lm_dist$m1 <- gsub("AP1", "AP-1", lm_dist$m1)
lm_dist$m2 <- gsub("AP1", "AP-1", lm_dist$m2)


### load in silico predictions
df_summary_main <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/DeepSTARR_motif_pair_distance_synth_sequences_df_summary_main.rds"))

candidate_motifs <- c(GATA="AGATAAGA",
                      "AP-1"="ATGACTCAT",
                      twist="ACATCTGT",
                      Trl="AGAGAGA",
                      ETS="ACCGGAAG",
                      SREBP="ATCACGCCA",
                      Dref="TATCGATA",
                      Ohler1="GTGTGACC",
                      Ohler6="AAAATACCA")

### plots per motifA at the center
pdf("Fig5B.pdf", width = 7, height = 7)
for(motif1 in names(candidate_motifs)){
  
  if(motif1 %in% c("Dref", "Ohler1", "Ohler6")){type="hk"}else{type="dev"}
  
  m1_strand <- "+"
  m2_strand <- "+"
  
  # table to use
  tmp_main <- df_summary_main[df_summary_main$Motif1 %in% motif1 & df_summary_main$Motif1_strand %in% m1_strand & df_summary_main$Motif2_strand %in% m2_strand,]
  
  # enhancer type
  if(type=="dev"){
    tmp_main <- tmp_main[!tmp_main$Motif2 %in% c("Dref", "Ohler1", "Ohler6"),]
    tmp_main$Backbone_motif2_act <- tmp_main$Backbone_motif2_act_dev
    tmp_main$Fold_change_BPNet_style <- tmp_main$Fold_change_BPNet_style_dev
    tmp_main$Residuals_FC <- tmp_main$Residuals_FC_dev
    tmp_main$Backbone_motif1_motif2_act <- tmp_main$Backbone_motif1_motif2_act_dev
  }
  if(type=="hk"){
    tmp_main <- tmp_main[tmp_main$Motif2 %in% c("Dref", "Ohler1", "Ohler6", "GGGCT"),]
    tmp_main$Backbone_motif2_act <- tmp_main$Backbone_motif2_act_hk
    tmp_main$Fold_change_BPNet_style <- tmp_main$Fold_change_BPNet_style_hk
    tmp_main$Residuals_FC <- tmp_main$Residuals_FC_hk
    tmp_main$Backbone_motif1_motif2_act <- tmp_main$Backbone_motif1_motif2_act_hk
  }
  
  # final summarised table
  tmp <- tmp_main %>%
    group_by(Distance_Motif2_Motif1, Motif2) %>%
    summarise(Backbone_motif2_act=median(Backbone_motif2_act),
              Fold_change_BPNet_style=median(Fold_change_BPNet_style),
              Residuals_FC=median(Residuals_FC),
              Backbone_motif1_motif2_act=median(Backbone_motif1_motif2_act))
  tmp$side <- factor(ifelse(tmp$Distance_Motif2_Motif1>0, "down", "up"), levels = c("up", "down"))
  
  # select pair motif
  if(type=="dev") m2_list <- c("GATA", "AP-1", "SREBP", "twist","ETS", "Trl")
  if(type=="hk") m2_list <- c("Dref", "Ohler1", "Ohler6")
  
  for(motif2 in m2_list){
    tmp2 <- tmp[tmp$Motif2 %in% motif2,]
    
    # Interaction residuals FC
    g_smooth_res <- ggplot(tmp2, aes(abs(Distance_Motif2_Motif1), Residuals_FC, col=Motif2, group=Motif2)) +
      geom_point(alpha=0.3)+
      geom_smooth(span = 0.1, se=F)+
      ggtitle(paste0(motif1, " / ", motif2)) +
      theme(plot.title = element_text(hjust=0.5)) +
      scale_y_continuous("DeepSTARR predicted coperativity") +
      scale_x_continuous("Motif pair distance (bp)", breaks = seq(-200,200,25), limits = c(0,125)) +
      geom_hline(yintercept = 1, linetype="dashed", col="grey60") +
      scale_color_manual("", values=high_col_list, drop=T) +
      guides(col=F) +
      theme(axis.title.x = element_blank(),
            axis.text = element_text(colour="black", size=12))
    
    ### add barplot of co-occurency
    fisher_results_tmp <- fisher_results[fisher_results$m1==motif1 & fisher_results$m2==motif2 & fisher_results$Distance %in% c("0-25", "25-50", "50-75", "75-100", "100-125"),]
    
    g_occurence <- ggplot(fisher_results_tmp, aes(Distance, log2(fisher_OR), fill=log2(fisher_OR))) +
      geom_bar(stat="identity") +
      geom_hline(yintercept = 0, linetype="dashed", col="grey60") +
      scale_fill_gradient2("Odds ratio", low = RColorBrewer::brewer.pal(7, "PRGn")[1], mid = "grey70", high = RColorBrewer::brewer.pal(7, "PRGn")[7], midpoint = 0) +
      # scale_fill_brewer("Odds ratio", palette="PRGn") +
      guides(fill=F) +
      ylab("Occurrence\n[log2 OR]") +
      xlab("Motif pair distance (bp)") +
      geom_text(aes(y=0, label=fisher_FDR_label), color="black", size=10) +
      theme(axis.text = element_text(colour="black", size=12),
            axis.title.x = element_text(colour="black", size=15))
    
    ### add barplot or co-occurency
    lm_dist_tmp <- lm_dist[lm_dist$m1==motif1 & lm_dist$m2==motif2 & lm_dist$cutoff %in% c("0-25", "25-50", "50-75", "75-100", "100-125"),]
    
    g_association <- ggplot(lm_dist_tmp, aes(cutoff, Estimate_lower_than, fill=Estimate_lower_than)) +
      geom_bar(stat="identity") +
      geom_hline(yintercept = 0, linetype="dashed", col="grey60") +
      scale_fill_gradient2("lm coef", low = "#2166AC", mid = "grey70", high = "#B2182B", midpoint = 0) +
      guides(fill=F) +
      scale_y_continuous("Enh activity\n[lm coef]") +
      xlab("Motif pair distance (bp)") +
      geom_text(aes(y=0, label=FDR_label), color="black", size=10) +
      theme(axis.title.x = element_blank(),
            axis.text = element_text(colour="black", size=12))
    
    print(g_smooth_res + g_association + g_occurence + plot_layout(ncol=1, heights = c(2.5,1,1)))
    
  }
  
  print(motif1)
}

dev.off()


########
# Fig 5C
########

### plots per motifA at the center
for(motif1 in names(candidate_motifs)){
  
  func="median"
  if(motif1 %in% c("Dref", "Ohler1", "Ohler6")){type="hk"}else{type="dev"}
  
  pdf(paste0("Fig5C_", as.character(motif1), ".pdf"), width = 7, height = 5.5)
  
  for(m1_strand in c("+", "-")){
    for(m2_strand in c("+", "-")){
      
      # table to use
      tmp_main <- df_summary_main[df_summary_main$Motif1 %in% motif1 & df_summary_main$Motif1_strand %in% m1_strand & df_summary_main$Motif2_strand %in% m2_strand,]
      
      # enhancer type
      if(type=="dev"){
        tmp_main$Backbone_motif2_act <- tmp_main$Backbone_motif2_act_dev
        tmp_main$Fold_change_BPNet_style <- tmp_main$Fold_change_BPNet_style_dev
        tmp_main$Residuals_FC <- tmp_main$Residuals_FC_dev
        tmp_main$Backbone_motif1_motif2_act <- tmp_main$Backbone_motif1_motif2_act_dev
      }
      if(type=="hk"){
        tmp_main$Backbone_motif2_act <- tmp_main$Backbone_motif2_act_hk
        tmp_main$Fold_change_BPNet_style <- tmp_main$Fold_change_BPNet_style_hk
        tmp_main$Residuals_FC <- tmp_main$Residuals_FC_hk
        tmp_main$Backbone_motif1_motif2_act <- tmp_main$Backbone_motif1_motif2_act_hk
      }
      
      # final summarised table
      tmp <- tmp_main %>%
        group_by(Distance_Motif2_Motif1, Motif2) %>%
        summarise(Backbone_motif2_act=median(Backbone_motif2_act),
                  Fold_change_BPNet_style=median(Fold_change_BPNet_style),
                  Residuals_FC=median(Residuals_FC),
                  Backbone_motif1_motif2_act=median(Backbone_motif1_motif2_act))
      tmp$side <- factor(ifelse(tmp$Distance_Motif2_Motif1>0, "down", "up"), levels = c("up", "down"))
      
      # Interaction residuals FC
      g_smooth_res <- ggplot(tmp, aes(abs(Distance_Motif2_Motif1), Residuals_FC, col=Motif2, group=Motif2)) +
        geom_point(alpha=0.3)+
        geom_smooth(span = 0.1, se=F)+
        ggtitle(paste0(motif1, m1_strand, "/motif2", m2_strand, " : ", func, " across ", length(unique(df_summary_main$Backbone)), " backbones")) +
        theme(plot.title = element_text(hjust=0.5)) +
        scale_y_continuous("DeepSTARR predicted coperativity") +
        scale_x_continuous("Motif pair distance (bp)", breaks = c(-10,10,seq(-200,200,25))) +
        geom_hline(yintercept = 1, linetype="dashed", col="grey60") +
        scale_color_manual("", values=high_col_list, drop=T)
      
      print(g_smooth_res)
      
    }
  }
  
  dev.off()
  
  print(motif1)
}


########
# Fig 5D,E
########

## Get confident motif positions

library(motifmatchr)
library(TFBSTools)
library(seqinr)

# Motifs
load(url("https://data.starklab.org/almeida/Motif_clustering/TF_clusters_PWMs.RData"))
View(TF_clusters_PWMs$metadata)
load(url("https://data.starklab.org/almeida/Drosophila_enhancers_motif_enrichment/Motif_enrichment_all.RData"))
View(Results_list_all$dev_vs_ctrl)
View(Results_list_all$hk_vs_ctrl)

TF_motifs <- list(srp=data.frame(Motif="GATA", ID="flyfactorsurvey__srp_SANGER_5_FBgn0003507", core="GATAA"),
                  kay_Jra=data.frame(Motif="AP1", ID="flyfactorsurvey__kay_Jra_SANGER_5_FBgn0001291", core="TGA.TCA"),
                  twist=data.frame(Motif="twist", ID="flyfactorsurvey__twi_da_SANGER_5_FBgn0000413", core="CAGATG"),
                  Trl=data.frame(Motif="Trl", ID="flyfactorsurvey__Trl_FlyReg_FBgn0013263", core="GAGA"),
                  ETS=data.frame(Motif="ETS", ID="flyfactorsurvey__Ets97D_SANGER_10_FBgn0004510", core="CCGGAA"),
                  SREBP=data.frame(Motif="SREBP", ID="flyfactorsurvey__HLH106_SANGER_10_FBgn0015234", core="TCACGCGA"),
                  Dref=data.frame(Motif="Dref", ID="homer__AVYTATCGATAD_DREF", core="TATCGATA"),
                  Ohler1=data.frame(Motif="Ohler1", ID="homer__MYGGTCACACTG_Unknown1", core="GGTCACACT"),
                  Ohler6=data.frame(Motif="Ohler6", ID="homer__AAAAATACCRMA_Unknown4", core="AAATACCA"))
TF_motifs <- do.call(rbind, TF_motifs)

# enhancers
twist_enhancers <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_mutation_individual_instances.rds"))
Final_enhancers_selected <- GRanges(paste0(sapply(strsplit(twist_enhancers$Sequence_ID,"_"), `[`, 1),
                                           ":", sapply(strsplit(twist_enhancers$Sequence_ID,"_"), `[`, 2),
                                           "-", sapply(strsplit(twist_enhancers$Sequence_ID,"_"), `[`, 3)),
                                    strand = twist_enhancers$Strand,
                                    Sequence_ID=twist_enhancers$Sequence_ID,
                                    Sequence=twist_enhancers$Sequence,
                                    enhancer_group=twist_enhancers$enhancer_group,
                                    seqinfo = Dmelanogaster@seqinfo)
Final_enhancers_selected <- unique(Final_enhancers_selected)

# motif positions in oligo
motif_ix_pos <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[name(TF_clusters_PWMs$All_pwms_log_odds) %in% TF_motifs$ID],
                            Final_enhancers_selected$Sequence,
                            genome = "BSgenome.Dmelanogaster.UCSC.dm3", p.cutoff = 5e-4, bg="genome", out = "positions")
names(motif_ix_pos) <- name(TF_clusters_PWMs$All_pwms_log_odds[name(TF_clusters_PWMs$All_pwms_log_odds) %in% TF_motifs$ID])

motif_ix_pos2 <- lapply(motif_ix_pos, function(motif_x){
  names(motif_x) <- Final_enhancers_selected$Sequence_ID
  motif_x <- motif_x[sapply(motif_x, length)>0] # remove sequences without motif
  # join all IRanges in same GRanges
  motif_x <- GRanges(names(unlist(motif_x)), #Use sequence ID as seqnames
                     IRanges(start(unlist(motif_x)), end(unlist(motif_x))),
                     strand = mcols(unlist(motif_x))$strand,
                     score=mcols(unlist(motif_x))$score)
  ### reduce - keep strand information
  # motif_x <- my_reduce(motif_x, min.frac.ov=0.7)
  ## instead, use plyranges::reduce_ranges to get motif score back - but here there is no min overlap
  # motif_x <- reduce_ranges(motif_x, max_score = max(score), sum_score = sum(score), Number=n()) # bad because one bp overlapping is too stringent, use the adapted version below
  motif_x <- my_reduce_with_score(motif_x, min.frac.ov=0.5)
  # add wt sequence
  motif_x$enh_strand <- twist_enhancers$Strand[match(as.character(motif_x@seqnames), twist_enhancers$Sequence_ID)]
  motif_x$seq <- substr(as.character(getSeq(Dmelanogaster, GRanges(paste0(sapply(strsplit(as.character(motif_x@seqnames),"_"), `[`, 1),
                                                                          ":", sapply(strsplit(as.character(motif_x@seqnames),"_"), `[`, 2),
                                                                          "-", sapply(strsplit(as.character(motif_x@seqnames),"_"), `[`, 3)),
                                                                   strand = motif_x$enh_strand))),
                        motif_x@ranges@start,
                        motif_x@ranges@start+motif_x@ranges@width-1)
  return(motif_x)
}) 
lapply(motif_ix_pos2, function(i) table(width(i)))

motif_ix_pos3 <- lapply(names(motif_ix_pos2), function(i){
  x <- data.frame(motif_ix_pos2[[i]], stringsAsFactors = F)
  names(x)[1] <- "Sequence"
  x$Sequence <- as.character(x$Sequence)
  x$strand <- as.character(x$strand)
  x$Motif <- TF_motifs$Motif[TF_motifs$ID %in% i]
  return(x)
})
names(motif_ix_pos3) <- names(motif_ix_pos2)

sapply(motif_ix_pos2, length)
sapply(motif_ix_pos3, nrow)

saveRDS(do.call(rbind, motif_ix_pos3), file = "Motif_oligo_positions_list.rds")



### motif pairs - associations with distance

# mutation data
mutation_data_2 <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_mutation_individual_instances_extended_info.rds"))
# mutation_data_2 <- mutation_data_2[-grep("all_mutant",mutation_data_2$Oligo_ID),]
mutation_data_2 <- mutation_data_2[complete.cases(mutation_data_2$enhancer_group),]

# confident TF motif positions
Motif_oligo_positions <- readRDS(url("Motif_oligo_positions_list.rds"))
Motif_oligo_positions$instance_center <- Motif_oligo_positions$start+(Motif_oligo_positions$end-Motif_oligo_positions$start)/2
table(Motif_oligo_positions$Motif)

TF_motifs <- list(srp=data.frame(Motif="GATA", ID="flyfactorsurvey__srp_SANGER_5_FBgn0003507", core="GATAA"),
                  kay_Jra=data.frame(Motif="AP1", ID="flyfactorsurvey__kay_Jra_SANGER_5_FBgn0001291", core="TGA.TCA"),
                  twist=data.frame(Motif="twist", ID="flyfactorsurvey__twi_da_SANGER_5_FBgn0000413", core="CAGATG"),
                  Trl=data.frame(Motif="Trl", ID="flyfactorsurvey__Trl_FlyReg_FBgn0013263", core="GAGA"),
                  ETS=data.frame(Motif="ETS", ID="flyfactorsurvey__Ets97D_SANGER_10_FBgn0004510", core="CCGGAA"),
                  SREBP=data.frame(Motif="SREBP", ID="flyfactorsurvey__HLH106_SANGER_10_FBgn0015234", core="TCACGCGA"),
                  Dref=data.frame(Motif="Dref", ID="homer__AVYTATCGATAD_DREF", core="TATCGATA"),
                  Ohler1=data.frame(Motif="Ohler1", ID="homer__MYGGTCACACTG_Unknown1", core="GGTCACACT"),
                  Ohler6=data.frame(Motif="Ohler6", ID="homer__AAAAATACCRMA_Unknown4", core="AAATACCA"))
TF_motifs <- do.call(rbind, TF_motifs)

for(m1 in TF_motifs$Motif[c(1:4,7)]){
  
  pdf(paste0("Fig5D_5E_", m1, ".pdf"), width = 3, height = 4.5)
  
  if(m1=="GATA") m1_twist <- "GATAA"
  if(m1=="AP1") m1_twist <- "TGA.TCA"
  if(m1=="Trl") m1_twist <- "GAGA"
  if(m1=="twist") m1_twist <- "CA..TG"
  if(m1=="Dref") m1_twist <- "ATCGAT"
  
  tmp <- Motif_oligo_positions[Motif_oligo_positions$Motif %in% m1 & Motif_oligo_positions$Sequence %in% mutation_data_2$Sequence_ID[mutation_data_2$Motif_mutated %in% m1_twist],]
  
  # mutations
  class <- ifelse(m1=="Dref", "hk", "dev")
  mutation_data_2_tmp <- mutation_data_2[mutation_data_2$Motif_mutated == m1_twist & mutation_data_2$enhancer_group %in% class,]
  
  tmp$log2FC <- sapply(1:nrow(tmp), function(i){
    if(m1!="Dref") out <- mutation_data_2_tmp$dev_log2FC_wt_mut[subjectHits(findOverlaps(GRanges(tmp$Sequence[i], IRanges(tmp$start[i], tmp$end[i])),
                                                                                         GRanges(mutation_data_2_tmp$Sequence_ID, IRanges(mutation_data_2_tmp$instance_start, mutation_data_2_tmp$instance_end)),
                                                                                         minoverlap = unique(nchar(mutation_data_2_tmp$wt_instance))))]
    if(m1=="Dref") out <- mutation_data_2_tmp$hk_log2FC_wt_mut[subjectHits(findOverlaps(GRanges(tmp$Sequence[i], IRanges(tmp$start[i], tmp$end[i])),
                                                                                        GRanges(mutation_data_2_tmp$Sequence_ID, IRanges(mutation_data_2_tmp$instance_start, mutation_data_2_tmp$instance_end)),
                                                                                        minoverlap = unique(nchar(mutation_data_2_tmp$wt_instance))))]
    if(length(out)>0){return(mean(out))}else{return(NA)}
  })
  tmp <- tmp[complete.cases(tmp$log2FC),]
  
  # get number of instances and their distance to each partner motif
  add_info <- lapply(TF_motifs$Motif, function(m2){
    Number <- sapply(1:nrow(tmp), function(s){
      m2_df <- Motif_oligo_positions[Motif_oligo_positions$Sequence %in% tmp$Sequence[s] & Motif_oligo_positions$Motif %in% m2,]
      if(nrow(m2_df)>0){return(nrow(m2_df))}else{return(0)}
    })
    Score <- sapply(1:nrow(tmp), function(s){
      m2_df <- Motif_oligo_positions[Motif_oligo_positions$Sequence %in% tmp$Sequence[s] & Motif_oligo_positions$Motif %in% m2,]
      if(nrow(m2_df)>0){return(sum(m2_df$max_score))}else{return(0)}
    })
    Distance <- sapply(1:nrow(tmp), function(s){
      m2_df <- Motif_oligo_positions[Motif_oligo_positions$Sequence %in% tmp$Sequence[s] & Motif_oligo_positions$Motif %in% m2 & !rownames(Motif_oligo_positions) %in% rownames(tmp)[s],]
      if(nrow(m2_df)>0){
        dist <- m2_df$instance_center-tmp$instance_center[s]
        return(dist[abs(dist)==min(abs(dist))][1])
      }else{return(NA)}
    })
    # double check
    out <- data.frame(Number, Score, Distance)
    names(out) <- paste0(m2, c("_n", "_score", "_dist"))
    return(out)
  })
  
  tmp2 <- cbind(tmp, do.call(cbind, add_info))
  
  # plots
  # select pair motif
  if(class=="dev") m2_list <- c("GATA", "AP1", "SREBP", "twist","ETS", "Trl")
  if(class=="hk") m2_list <- c("Dref", "Ohler1", "Ohler6")
  
  for(m2 in m2_list){
    
    gg_df <- tmp2[complete.cases(tmp2[,paste0(m2, "_dist")]),] # restricting by complete distances means it has at least 2 when is homotypic pair, and at least one of each when heterotypic. Should I limit homotypic to 2 instances
    
    # remove overlapping instances (distance between centers > motif length)
    gg_df <- gg_df[abs(gg_df[,paste0(m2, "_dist")]) > median(gg_df$width),]
    
    gg_df$class <- NA
    gg_df$class[abs(gg_df[,paste0(m2, "_dist")])<25] <- "<25"
    gg_df$class[abs(gg_df[,paste0(m2, "_dist")])>50] <- ">50"
    table(gg_df$class)
    
    g2 <- ggplot(gg_df[complete.cases(gg_df$class),], aes(class, log2FC, colour=class)) +
      geom_boxplot(fill=NA, size=1) +
      scale_colour_manual(values=c("<25"="#2166AC", ">50"="#B2182B")) +
      guides(col=F) +
      scale_y_continuous(paste0("log2 FC enhancer activity [",m1,"]"), breaks = seq(-10,10,2)) +
      scale_x_discrete(paste0(m1, "-", m2," distance (bp)"),
                       labels=c(paste0("<25\n(n=", length(which(gg_df$class=="<25")),")"),
                                paste0(">50\n(n=", length(which(gg_df$class==">50")),")"))) +
      geom_hline(yintercept = 0, linetype="dashed", col="grey60")
    
    wilcox.test(gg_df$log2FC~gg_df$class)
    
    print(g2)
    
  }
  
  dev.off()
  print(m1)
  
}

