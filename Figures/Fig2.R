
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(GenomicRanges)
library(ggplot2)
theme_set(theme_light() + theme(axis.text = element_text(colour = "black")))
library(ggpointdensity)

########
# Fig 2A
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
pdf("Fig2A.pdf", height = 2.1, width = 13)
# dev
i="chr3L_3310914_3311162_-_wt_dCP"
my.logo(twist_contr_scores$dev[[grep(i, names(twist_contr_scores$dev), fixed = T)]]) + ggtitle(paste0(i, " - dev scores"))

# hk
i="chrX_4580794_4581042_-_wt_hkCP"
my.logo(twist_contr_scores$hk[[grep(i, names(twist_contr_scores$hk), fixed = T)]]) + ggtitle(paste0(i, " - dev scores"))

dev.off()


########
# Fig 2C
########

df <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_mutation_all_instances_results.txt")
# only mutant version 1
df <- df[df$Mutant_version %in% "s1",]

# per enhancer type
boxplot_list <- list()
for(CP in c("dev", "hk")){
  
  # only strong and specific dev or hk enhancers
  if(CP=="dev") df_CP <- df[complete.cases(df$dev_mut) & df$enhancer_group %in% "dev",]
  if(CP=="hk") df_CP <- df[complete.cases(df$hk_mut) & df$enhancer_group %in% "hk",]
  
  # rm GATA (in hk because it is similar to Dref) or GATA_no_DRE (in dev)
  if(CP=="dev"){
    df_CP <- df_CP[!df_CP$Motif_mutated %in% c("GATA (no Dref)"),]
  }
  if(CP=="hk"){
    df_CP <- df_CP[!df_CP$Motif_mutated %in% c("GATA"),]
  }
  
  ### rm motifs with only one enhancer
  df_CP <- df_CP[!df_CP$Motif_mutated %in% levels(df_CP$Motif_mutated)[table(df_CP$Motif_mutated)<2],]
  df_CP$Motif_mutated <- droplevels(df_CP$Motif_mutated)
  
  df_CP$Motif_mutated2 <- as.character(df_CP$Motif_mutated)
  df_CP$Motif_mutated2[grep("ctrl", df_CP$Motif_mutated)] <- "3 controls"
  df_CP$Motif_mutated2 <- factor(df_CP$Motif_mutated2,
                                 levels = c("3 controls", "AP-1", "GATA", "GATA (no Dref)", "SREBP", "CREB", "twist", "ETS", "STAT", "Trl", "Dref", "Ohler1", "Ohler7", "Ohler6"))
  
  df_CP$Motif_mutated3 <- as.character(df_CP$Motif_mutated)
  df_CP$Motif_mutated3[-grep("ctrl", df_CP$Motif_mutated)] <- "Motif"
  df_CP$Motif_mutated3 <- factor(df_CP$Motif_mutated3,
                                 levels = c("ctrl_TAGG", "ctrl_CCTTA", "ctrl_GGGCT", "Motif"))
  
  ### boxplot all motifs
  
  # x labels counts per boxplot
  xlabels <- sapply(as.character(levels(df_CP$Motif_mutated2)), function(x){
    paste0(gsub("_", " ", x)," (", table(df_CP$Motif_mutated2)[x], ")")
  })
  
  motif_colours <- sapply(names(xlabels), function(i){
    if(length(grep("3 controls", i))>0) return("grey60")
    if(i %in% c("GATA", "GATA (no Dref)", "AP-1", "GAGA", "Trl", "twist", "SREBP", "ETS", "STAT", "CREB")) return("orangered")
    if(i %in% c("CAGCTG", "Dref", "Ohler1", "Ohler6", "Ohler7")) return("dodgerblue")
  })
  
  # if less than 5 points, plot points instead of boxplot
  df_CP$FC <- df_CP[,paste0(CP, "_log2FC_wt_mut")]
  boxplot <- ggplot(df_CP, aes(Motif_mutated2, FC, fill=Motif_mutated2, alpha=Motif_mutated3)) +
    geom_point(size=-1) +
    geom_boxplot(data=df_CP[!df_CP$Motif_mutated2 %in% levels(df_CP$Motif_mutated2)[table(df_CP$Motif_mutated2)<=5],], aes(Motif_mutated2, FC, fill=Motif_mutated2, alpha=Motif_mutated3),
                 outlier.size = 0.6) +
    geom_boxplot(data=df_CP[df_CP$Motif_mutated2 %in% levels(df_CP$Motif_mutated2)[table(df_CP$Motif_mutated2)<=5],], aes(Motif_mutated2, FC, fill=Motif_mutated2, alpha=Motif_mutated3),
                 outlier.size = 0.6, alpha=0.2) +
    geom_point(data=df_CP[df_CP$Motif_mutated2 %in% levels(df_CP$Motif_mutated2)[table(df_CP$Motif_mutated2)<=5],], aes(Motif_mutated2, FC, fill=Motif_mutated2, alpha=Motif_mutated3),
               size=2.2, shape=21, position = position_jitter(width = 0.15)) +
    ylab(paste0("log2 FC ", CP, " enhancer activity")) +
    scale_x_discrete("Mutated motifs", breaks= names(xlabels), labels = xlabels) +
    geom_hline(yintercept = 0, linetype="dashed", col="grey40") +
    scale_fill_manual(values=motif_colours) +
    scale_alpha_manual(values=rep(1,5)) +
    guides(fill=F, alpha=F) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(colour="black", size=15, angle=45, hjust=1),
          axis.text.y = element_text(colour="black", size=14),
          axis.title = element_text(colour="black", size=16)) +
    scale_y_continuous(breaks = seq(-10,20,2), limits = c(-8, 3))
  
  boxplot_list[[CP]] <- boxplot
  
}

pdf("Fig2C.pdf", width = 6.5, height = 5.5)
print(boxplot_list$dev)
print(boxplot_list$hk)
dev.off()


########
# Fig 2D
########

library(ggrepel)
library(dplyr)

high_col_list_tmp <- c(GATA=rgb(116,159,242, maxColorValue = 255),
                       "AP-1"=rgb(177,53,115, maxColorValue = 255),
                       AP1=rgb(177,53,115, maxColorValue = 255),
                       Trl=rgb(212,147,91, maxColorValue = 255),
                       twist=rgb(121,170,109, maxColorValue = 255),
                       ETS=rgb(255,102,102, maxColorValue = 255),
                       SREBP=rgb(102,102,0, maxColorValue = 255),
                       CREB3="#CCCC00",
                       MAF="#E41A1C",
                       Dref="#51A9FF",
                       Ohler1="#1E65AB",
                       Ohler5="#4C0099",
                       "Ebox/Ohler5"="#4C0099",
                       Ohler6="#24ACAC",
                       Ohler7="#B266FF")

df <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/DeepSTARR_motif_imp_and_motif_enrichment.txt")

# highlight specific groups
df$motif_group2 <- NA
df$motif_group2[grep("^AP1", df$Motif_cluster_name, ignore.case = T)] <- "AP-1"
df$motif_group2[grep("^CREB3", df$Motif_cluster_name, ignore.case = T)] <- "CREB3"
df$motif_group2[grep("^GATA", df$Motif_cluster_name, ignore.case = T)] <- "GATA"
df$motif_group2[grep("^Trl|GAGA-repeat", df$Motif_cluster_name, ignore.case = T)] <- "Trl"
df$motif_group2[grep("ETS/4|ETS/5|ETS/1|ETS/7|ETS/2|^ETS$", df$Motif_cluster_name, ignore.case = T)] <- "ETS"
df$motif_group2[grep("SREBP", df$Motif_cluster_name, ignore.case = T)] <- "SREBP"
df$motif_group2[grep("^twi", df$motif_description2, ignore.case = T)] <- "twist"
df$motif_group2[grep("^Ebox/CATCTG$", df$motif_description2, ignore.case = T)] <- "twist"

df$motif_group2[grep("DRE", df$Motif_cluster_name, ignore.case = T)] <- "Dref"
df$motif_group2[grep("Ohler1", df$Motif_cluster_name, ignore.case = T)] <- "Ohler1"
df$motif_group2[grep("Ohler6", df$Motif_cluster_name, ignore.case = T)] <- "Ohler6"
df$motif_group2[grep("^Ohler7", df$Motif_cluster_name, ignore.case = T)] <- "Ohler7"

df$motif_group2[grep("dev_new", df$Motif)][!complete.cases(df$motif_group2[grep("dev_new", df$Motif)])] <- "Others dev"
df$motif_group2[grep("hk_new", df$Motif)][!complete.cases(df$motif_group2[grep("hk_new", df$Motif)])] <- "Others hk"
df$motif_group2[grep("DRE_Ohler7", df$motif_description2)] <- "Others hk"

theme_set(theme_classic(base_size=14) + theme(axis.text = element_text(colour = "black")))

pdf("Fig2D.pdf", width = 7.5, height = 5.5)
for(class in c("dev", "hk")){
  
  if(class=="dev"){
    df$X <- df$Enrichment_dev_enhancers_log2OR
    df$Y <- df$log2FC_Dev
    tmp <- df[complete.cases(df$motif_group2) & !df$motif_group2 %in% c("Dref", "Ohler1", "Ohler6", "Ohler7", "Ebox/Ohler5", "Others hk"),]
  }else if(class=="hk"){
    df$X <- df$Enrichment_hk_enhancers_log2OR
    df$Y <- df$log2FC_Hk
    tmp <- df[complete.cases(df$motif_group2) & df$motif_group2 %in% c("Dref", "Ohler1", "Ohler6", "Ohler7", "Ebox/Ohler5", "Others hk"),]
  }
  
  # Cluster representatives (two best by motif enrichment and DeepSTARR prediction)
  tmp2 <- rbind(df %>%
                  group_by(Motif_cluster) %>%
                  top_n(1, X),
                df %>%
                  group_by(Motif_cluster) %>%
                  top_n(1, Y))
  tmp2 <- tmp2[!duplicated(tmp2$Motif),]
  
  gg_type <- ggplot(tmp2, aes(X, Y)) +
    geom_point(col="grey70", size=0.8) +
    scale_x_continuous("Motif enrichment (log2 odds ratio)", breaks = seq(-10,10,1)) +
    scale_y_continuous("DeepSTARR importance (log2 FC)", breaks = seq(-10,10,0.5)) +
    geom_hline(yintercept = 0, linetype="dashed", col="grey60") +
    geom_vline(xintercept = 0, linetype="dashed", col="grey60") +
    geom_point(data=tmp2[tmp2$motif_group2 %in% tmp$motif_group2,], aes(X, Y, col=motif_group2),
               size=2.2) +
    scale_colour_manual("Motif group", values = c(high_col_list_tmp, "Others dev"="black", "Others hk"="black")) +
    theme(axis.title=element_text(size=18),
          axis.text=element_text(size=16),
          legend.title=element_text(size=18),
          legend.text=element_text(size=16))
  
  print(gg_type)
  
}
dev.off()


