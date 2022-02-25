
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(GenomicRanges)
library(ggplot2)
theme_set(theme_light() + theme(axis.text = element_text(colour = "black")))
library(ggpointdensity)
library(stringdist)
library(cowplot)

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
# Fig 4A
########

# load sequences
Peaks <- data.frame(readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_dev_hk_neg_specific_enhancers_information.rds")), stringsAsFactors = F)
Peaks$Class <- as.character(Peaks$Class)
table(Peaks$Class, Peaks$enh_type)

# load contribution scores
contr_scores <- readRDS(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/DeepSTARR_contr_scores_dev_hk_neg_enhancers.rds"))

motif_interest_list_dev <- c("GATAA", "TGA.TCA", "CATCTG", "GAGAG")
motif_interest_list_hk <- c("ATCGAT")
motif_interest_list_ctrl <- c("GGGCT")
motif_list <- unique(c(motif_interest_list_dev, motif_interest_list_hk, motif_interest_list_ctrl))
names(motif_list) <- motif_list

label_list <- c("dev", "hk")
names(label_list) <- label_list

motif_scores_list <- lapply(c("flank5"=5, "flank10"=10, "flank50"=50), function(flank_size){
  
  return(lapply(as.list(label_list), function(class){
    
    class_scores <- contr_scores[[class]]
    
    if(class == "dev") motif_list2 <- c(GATAA="GATAA", "TGA.TCA"="TGA.TCA",  CATCTG="CATCTG", GAGAG="GAGAG", GGGCT="GGGCT")
    if(class == "hk") motif_list2 <- c(ATCGAT="ATCGAT",  GGGCT="GGGCT")
    
    motif_scores_list_class <- lapply(motif_list2, function(motif){
      
      pos_fw <- stringr::str_locate_all(Peaks$Sequence, motif)
      names(pos_fw) <- Peaks$ID
      pos_fw <- pos_fw[sapply(pos_fw, length)>0]
      pos_fw <- cbind(dplyr::bind_rows(lapply(pos_fw, as.data.frame), .id="Sequence_ID"), data.frame(strand="fw", stringsAsFactors = F))
      
      pos_rv <- stringr::str_locate_all(Peaks$Sequence, as.character(reverseComplement(DNAString(motif))))
      names(pos_rv) <- Peaks$ID
      pos_rv <- pos_rv[sapply(pos_rv, length)>0]
      pos_rv <- cbind(dplyr::bind_rows(lapply(pos_rv, as.data.frame), .id="Sequence_ID"), data.frame(strand="rv", stringsAsFactors = F))
      
      pos <- rbind(pos_fw, pos_rv)[,c(2:4,1)]
      # remove rev motif occurences if motif is palindromic
      # when motif is palindromic but alows N (e.g. CA..TG), I only keep the fw seqeunce anyway
      if(motif==as.character(reverseComplement(DNAString(motif)))) pos <- pos_fw[,c(2:4,1)]
      
      motif_scores <- do.call(rbind, lapply(1:nrow(pos), function(s){
        
        # contribution scores
        ID <- pos$Sequence_ID[s]
        seq <- Peaks$Sequence[match(ID, Peaks$ID)]
        tmp <- class_scores[[grep(ID, names(class_scores), fixed = T)]]
        
        x <- pos[s,]
        x[,1] <- x[,1]-flank_size # extend flanks
        x[,2] <- x[,2]+flank_size # extend flanks
        
        scores <- sapply(x[,1]:x[,2], function(n){
          if(n<1 | n>249){return(NA)}else{return(sum(tmp[,n]))}
        })
        if(x[,3] == "rv") rev(scores)
        names(scores) <- c(-flank_size:-1, sapply(1:nchar(motif), function(a) substr(motif, a,a)), 1:flank_size)
        
        inst_seq <- substr(seq, as.numeric(x[,1]), as.numeric(x[,2]))
        if(x[,3] == "rv") inst_seq <- as.character(reverseComplement(DNAString(inst_seq)))
        
        # mean of defined nucleotides (do not count the ".")
        if(length(grep("\\.", motif))==0) Avg <- mean(colSums(tmp[,pos[s,1]:pos[s,2]]))
        if(length(grep("\\.", motif))>0) Avg <- mean(colSums(tmp[,pos[s,1]:pos[s,2]][,-grep("\\.", sapply(1:nchar(motif), function(a) substr(motif, a,a)))]))
        
        c(as.character(pos[s,]),
          scores,
          Sequence=inst_seq,
          Average_core=Avg,
          range_core_nucleotides=range(colSums(tmp[,pos[s,1]:pos[s,2]]))[2]-range(colSums(tmp[,pos[s,1]:pos[s,2]]))[1])
      }))
      
      motif_scores_all <- as.data.frame(motif_scores)
      names(motif_scores_all)[1:4] <- c("start", "end", "strand", "Sequence_ID")
      motif_scores_all[,c(1,2,5:(ncol(motif_scores_all)-3), (ncol(motif_scores_all)-1):(ncol(motif_scores_all)))] <- apply(motif_scores_all[,c(1,2,5:(ncol(motif_scores_all)-3), (ncol(motif_scores_all)-1):(ncol(motif_scores_all)))], 2, function(x) as.numeric(as.character(x)))
      motif_scores_all$class <- Peaks$Class[match(motif_scores_all$Sequence_ID, Peaks$ID)]
      
      return(motif_scores_all)
    })
    
    return(motif_scores_list_class)
    
  }))
})

save(motif_scores_list, file="Motif_word_and_flanks_importances.RData")



### plot importance of flanking nucleotides
col_scheme <- c(A='#109648', "A.1"='#109648', "A.2"='#109648', "A.3"='#109648',
                C='#255C99', "C.1"='#255C99', "C.2"='#255C99',
                G='#F7B32B', "G.1"='#F7B32B', "G.2"='#F7B32B', "G.3"='#F7B32B',
                "T"='#D62839', "T.1"='#D62839', "T.2"='#D62839',
                "."="grey80")

label_list <- c("dev", "hk")
names(label_list) <- label_list

load(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophila_motif_word_and_flanks_importances.RData"))
flank_size <- "flank50"

pdf("Fig4A.pdf", height = 3.5, width = 8.5)
for(class in label_list){
  
  if(class == "dev") motif_list2 <- c("GATAA", "TGA.TCA",  "CATCTG",   "GAGAG",  "GGGCT")
  if(class == "hk") motif_list2 <- c("ATCGAT",  "GGGCT")
  
  for(motif in motif_list2){
    
    test <- motif_scores_list[[flank_size]][[class]][[motif]]
    # enhancer type
    test <- test[test$class %in% class,]
    
    # keep the ones with at least 10bp flanks
    test <- test[nchar(as.character(test$Sequence)) >= 20+nchar(motif),]
    head(sort(table(substr(test$Sequence, 51,50+nchar(motif))), decreasing = T))
    
    plot_list <- lapply(c(0.9), function(top){
      
      # top instances
      df_top <- reshape2::melt(test[test$Average_core >=quantile(test$Average_core, top), c(5,15,25,35,45:(5+49+nchar(motif)+10),74+nchar(motif), 84+nchar(motif), 94+nchar(motif), 104+nchar(motif))])
      df_tmp <- data.frame(variable=c("-45", "-35", "-25", "-15", "15", "25", "35", "45"),
                           value=NA, stringsAsFactors = F)
      df_top2 <- rbind(df_top, df_tmp)
      df_top2$variable <- factor(df_top2$variable, levels = c(-50, -45, -40, -35, -30, -25, -20, -15,
                                                              as.character(unique(df_top2$variable)[5:(24+nchar(motif))]),
                                                              paste0(seq(15,50,5))))
      df_top_for_line <- df_top %>% group_by(variable) %>% summarise(median = median(value, na.rm=T))
      
      # bottom instances
      df_bottom <- reshape2::melt(test[test$Average_core <=quantile(test$Average_core, 1-top), c(5,15,25,35,45:(5+49+nchar(motif)+10),74+nchar(motif), 84+nchar(motif), 94+nchar(motif), 104+nchar(motif))])
      df_bottom2 <- rbind(df_bottom, df_tmp)
      df_bottom2$variable <- factor(df_bottom2$variable, levels = c(-50, -45, -40, -35, -30, -25, -20, -15,
                                                                    as.character(unique(df_bottom2$variable)[5:(24+nchar(motif))]),
                                                                    paste0(seq(15,50,5))))
      df_bottom_for_line <- df_bottom %>% group_by(variable) %>% summarise(median = median(value, na.rm=T))
      
      # wilcoxon per position
      df_top_all <- reshape2::melt(test[test$Average_core >=quantile(test$Average_core, top), 5:(5+49+nchar(motif)+50)])
      df_top_all_for_line <- df_top_all %>% group_by(variable) %>% summarise(median = median(value, na.rm=T))
      df_bottom_all <- reshape2::melt(test[test$Average_core <=quantile(test$Average_core, 1-top), 5:(5+49+nchar(motif)+50)])
      df_bottom_all_for_line <- df_bottom_all %>% group_by(variable) %>% summarise(median = median(value, na.rm=T))
      pv <- sapply(unique(df_top_all$variable), function(i) wilcox.test(df_top_all$value[df_top_all$variable==i], df_bottom_all$value[df_bottom_all$variable==i])$p.value )
      pv <- ifelse(pv < 0.001, "*", "")
      names(pv) <- unique(df_top_all$variable)
      ylabel <- sapply(unique(df_top_all$variable), function(i) quantile(df_top_all$value[df_top_all$variable==i], 0.75, na.rm=T) + 1.8*IQR(df_top_all$value[df_top_all$variable==i], na.rm=T))
      names(ylabel) <- unique(df_top_all$variable)
      
      # sum of importance in flanks and core
      delta <- df_top_all_for_line$median-df_bottom_all_for_line$median
      names(delta) <- df_top_all_for_line$variable
      delta <- delta[pv=="*"]
      flankL <- sum(delta[1:length(which(pv[1:(grep("-1$", names(pv)))]=="*"))])
      flankR <- sum(delta[(grep("^1$", names(pv[pv=="*"]))):length(delta)])
      core <- sum(delta)-flankL-flankR
      
      gg_line <- ggplot() +
        geom_boxplot(data=df_top2, aes(factor(variable), value, fill="Top 90th"), outlier.size = -1, alpha=0.2, size=0.3, position = position_nudge(x = -0.2, y = 0), width=0.4) +
        geom_path(data = df_top_for_line, aes(variable, median, group = 1), col=RColorBrewer::brewer.pal(3, "Dark2")[3], size=0.8, position = position_nudge(x = -0.2, y = 0)) +
        geom_boxplot(data=df_bottom2, aes(factor(variable), value, fill="Bottom 10th"), outlier.size = -1, alpha=0.2, size=0.3, position = position_nudge(x = 0.2, y = 0), width=0.4) +
        geom_path(data = df_bottom_for_line, aes(variable, median, group = 1), col=RColorBrewer::brewer.pal(3, "Dark2")[1], size=0.8, position = position_nudge(x = 0.2, y = 0)) +
        scale_x_discrete(labels=c(-50, "", -40, "", -30, "", -20, "", -10, "", -8, "", -6, "", -4, "", -2, "",
                                  strsplit(motif, "")[[1]],
                                  "", "+2", "", "+4", "", "+6", "", "+8", "", "+10", "", "+20", "", "+30", "", "+40", "", "+50")) +
        ylab(paste0("DeepSTARR nucleotide importance")) +
        ggtitle(paste0(motif, " instances (", class, ")")) +
        scale_fill_manual(paste0("Instances"), values = c("Top 90th"=RColorBrewer::brewer.pal(3, "Dark2")[3],
                                                          "Bottom 10th"=RColorBrewer::brewer.pal(3, "Dark2")[1])) +
        guides(fill = guide_legend(reverse = T)) +
        geom_text(data=data.frame(x=names(pv)[names(pv) %in% unique(df_top$variable)],
                                  y=ylabel[names(ylabel) %in% unique(df_top$variable)],
                                  pv=pv[names(pv) %in% unique(df_top$variable)]),
                  aes(x,y, label=pv), color="black", size=4, inherit.aes = F) +
        geom_vline(xintercept = c(18.5, 18.5+nchar(motif)), linetype="dashed", col="grey50") +
        annotate("text", label=paste0("core\n(imp=",round(core,2),")"), x=18.5+nchar(motif)/2, y=Inf, vjust=1.2, color="grey30", size=3.5, lineheight=0.8) +
        annotate("text", label=c(paste0("flank up\n(imp=",round(flankL,2),")"), paste0("flank down\n(imp=",round(flankR,2),")")),
                 x=c(14.5, 22.5+nchar(motif)), y=Inf, vjust=1.2, color="grey30", size=3.5, lineheight=0.8) +
        theme(panel.background = element_rect(fill="white",colour="black"),
              axis.ticks.y = element_line(colour = "black"),
              axis.ticks = element_line(colour = c(rep(c("black", "transparent"),4),
                                                   rep("black", nchar(motif)+20),
                                                   rep(c("transparent", "black"),4))),
              axis.text=element_text(size=10, colour="black"), axis.title=element_text(size=14, colour="black"),
              axis.title.x = element_blank(),
              legend.key=element_rect(fill=NA, colour=NA),
              legend.text = element_text(size=11, colour="black"),
              legend.title = element_text(size=12, colour="black"),
              plot.title = element_text(size=16, hjust = 0.5, colour="black"), plot.subtitle = element_text(size=14, hjust = 0.5))
      
      if(motif=="GATAA"){
        gg_line <- gg_line + coord_cartesian(ylim=c(-0.05, max(ylabel)+0.05))
      }
      if(motif=="TGA.TCA"){
        gg_line <- gg_line + coord_cartesian(ylim=c(-0.08, max(ylabel)+0.05))
      }
      if(motif=="CATCTG"){
        gg_line <- gg_line + coord_cartesian(ylim=c(-0.07, max(ylabel)+0.05))
      }
      if(motif=="GAGAG"){
        gg_line <- gg_line + coord_cartesian(ylim=c(-0.05, max(ylabel)+0.03))
      }
      if(motif=="ATCGAT"){
        gg_line <- gg_line + coord_cartesian(ylim=c(-0.1, max(ylabel)+0.1))
      }
      
      return(gg_line)
    })
    
    print(plot_grid(plotlist = plot_list, nrow = 1))
    
  }
}
dev.off()


########
# Fig 4B
########

motif_interest_list_dev <- c("GATAA", "TGA.TCA", "CATCTG", "GAGAG")
motif_interest_list_hk <- c("ATCGAT")
motif_interest_list_ctrl <- c("GGGCT")
motif_list <- unique(c(motif_interest_list_dev, motif_interest_list_hk, motif_interest_list_ctrl))
names(motif_list) <- motif_list

label_list <- c("dev", "hk")
names(label_list) <- label_list

# run in looping
motif_scores_list <- lapply(as.list(label_list), function(class){
  
  class_scores <- contr_scores[[class]]
  
  motif_scores_list_class <- lapply(as.list(motif_list), function(motif){
    
    motif_scores <- lapply(1:nrow(Peaks), function(s){
      
      ID <- Peaks$ID[s]
      seq <- Peaks$Sequence[s]
      pos <- stringr::str_locate_all(seq, c(motif, as.character(reverseComplement(DNAString(motif)))))
      names(pos) <- c("fw", "rv")
      pos <- pos[lapply(pos,length)>0] # remove empty elements
      
      # contribution scores
      tmp <- class_scores[[grep(ID, names(class_scores), fixed = T)]]
      
      if(length(pos)>0){
        
        # remove rev motif occurences if motif is palindromic
        # when motif is palindromic but alows N (e.g. CA..TG), I only keep the fw seqeunce anyway
        if(motif==as.character(reverseComplement(DNAString(motif)))) pos <- pos[-2]
        
        pos <- lapply(1:length(pos), function(i){
          x <- pos[[i]]
          x <- cbind(x, matrix(rep(NA, nrow(x)*(nchar(motif)+7)), nrow = nrow(x))) # add space for nucleotide scores and average
          x[,3] <- names(pos)[i]
          for(r in 1:nrow(x)){
            if(names(pos)[i] == "fw") x[r,4:(3+nchar(motif))] <- colSums(tmp[,x[r,1]:x[r,2]])
            if(names(pos)[i] == "rv") x[r,4:(3+nchar(motif))] <- rev(colSums(tmp[,x[r,1]:x[r,2]])) # if reverse strand, reverse values to have that in the same position as the motif word
            # save the flanks
            flank_size=5
            if(names(pos)[i] == "fw") x[r,(3+nchar(motif)+1)] <- substr(seq, as.numeric(x[r,1])-flank_size, as.numeric(x[r,2])+flank_size)
            if(names(pos)[i] == "rv") x[r,(3+nchar(motif)+1)] <- as.character(reverseComplement(DNAString(substr(seq, as.numeric(x[r,1])-flank_size, as.numeric(x[r,2])+flank_size))))
            
            # mean of defined nucleotides (do not count the ".")
            if(length(grep("\\.", sapply(1:nchar(motif), function(a) substr(motif, a,a))))==0) x[r,(3+nchar(motif)+2)] <- mean(colSums(tmp[,x[r,1]:x[r,2]]))
            if(length(grep("\\.", sapply(1:nchar(motif), function(a) substr(motif, a,a))))>0) x[r,(3+nchar(motif)+2)] <- mean(colSums(tmp[,x[r,1]:x[r,2]][,-grep("\\.", sapply(1:nchar(motif), function(a) substr(motif, a,a)))]))
            
            x[r,(3+nchar(motif)+3)] <- range(colSums(tmp[,x[r,1]:x[r,2]]))[2]-range(colSums(tmp[,x[r,1]:x[r,2]]))[1] # range_nucleotides
            x[r,(3+nchar(motif)+4)] <- ID # Sequence_ID
            x[r,(3+nchar(motif)+5)] <- names(class_scores)[grep(ID, names(class_scores), fixed = T)] # Sequence_ID_long
            x[r,(3+nchar(motif)+6)] <- Peaks$Class[s] # class
            colnames(x) <- c("start", "end", "strand", sapply(1:nchar(motif), function(a) substr(motif, a,a)), "Sequence", "Average", "range_nucleotides", "Sequence_ID", "Sequence_ID_long", "Sequence_class")
          }
          
          return(x)
        })
        
        return(do.call(rbind, pos))
      }
    })
    
    ## add how many motifs per sequence
    motif_scores <- motif_scores[lapply(motif_scores,length)>0] # remove empty elements
    if(length(motif_scores)>0){
      motif_scores <- lapply(1:length(motif_scores), function(i){
        x <- motif_scores[[i]]
        x <- cbind(x,
                   data.frame(Range_average=range(as.numeric(x[,"Average"]))[2]-range(as.numeric(x[,"Average"]))[1],
                              Motif_counts=nrow(x)))
        return(x)
      })
      motif_scores_all <- as.data.frame(do.call(rbind, motif_scores))
      motif_scores_all[,c(1,2,4:(ncol(motif_scores_all)-8), (ncol(motif_scores_all)-6):(ncol(motif_scores_all)-5))] <- apply(motif_scores_all[,c(1,2,4:(ncol(motif_scores_all)-8), (ncol(motif_scores_all)-6):(ncol(motif_scores_all)-5))], 2, function(x) as.numeric(as.character(x)))
      
      return(motif_scores_all)
    }
    print(paste0(class, " >> ", motif))
  })
  
  return(motif_scores_list_class)
  
})

save(motif_scores_list, file="Motif_word_importances.RData")


### plots

# load importance scores of motifs
load(url("https://data.starklab.org/almeida/DeepSTARR/Figures_data/Drosophla_motif_word_importances.RData"))

motif_interest_list2 <- c("GATAA", "TGA.TCA", "CATCTG", "GAGAG", "ATCGAT", "GGGCT")

col_scheme <- c('#109648', '#255C99', '#F7B32B', '#D62839')

label_list <- c("dev", "hk")
names(label_list) <- label_list

for(class in label_list){
  
  if(class == "dev") my_col <- "orangered"
  if(class == "hk") my_col <- "dodgerblue"
  
  pdf(paste0("Fig4B_", class, ".pdf"), height = 8, width = 8)
  for(motif in motif_interest_list2){
    
    test <- motif_scores_list[[class]][[motif]]
    # enhancer type
    test <- test[test$Sequence_class %in% class,]
    
    # keep the ones with same flank size
    test <- test[nchar(as.character(test$Sequence))==10+nchar(motif),]
    table(substr(test$Sequence, 6,5+nchar(motif)))
    
    df<- data.frame()
    for(i in c(1:(nchar(motif)+10))){
      df <- rbind(df, data.frame(Sequence=paste(test$Sequence_ID, test$start, test$end, sep="_"),
                                 Pos=i,
                                 Base=substr(test$Sequence, i,i),
                                 Av.motif.imp=test$Average))
    }
    
    df$Pos <- factor(df$Pos, levels = c(1:(nchar(motif)+10)),
                     labels=c(1:5,
                              paste0("M",1:nchar(motif)),
                              (nchar(motif)+6):(nchar(motif)+10)))
    
    # order motifs
    df <- df[order(df$Av.motif.imp),]
    df$Sequence <- factor(df$Sequence, levels = unique(df$Sequence))
    
    heatmap <- ggplot(df, aes(factor(Pos), Sequence)) +
      geom_tile(aes(fill = Base), alpha=0.8) +
      scale_fill_manual(values = col_scheme) +
      scale_y_discrete(paste0(nrow(test), " instances ordered by average importance"), expand=c(0,0)) +
      scale_x_discrete("Flanking nucleotides", labels=c(-5:-1,
                                                        strsplit(motif, "")[[1]],
                                                        paste0("+", 1:5))) +
      theme(panel.background = element_rect(fill="white",colour="white"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.text.x = element_text(size = 12, colour="black"),
            axis.text.y = element_blank(),
            axis.title=element_text(size=14, colour="black"),
            axis.title.x = element_blank(),
            plot.title = element_text(size=14, hjust = 0.5, colour="black"),
            legend.key=element_rect(fill=NA, colour=NA),
            legend.text = element_text(size=11, colour="black"),
            legend.title = element_text(size=12, colour="black"))
    
    # print(heatmap)
    
    # statistical test for each nucleotide per position
    ids <- as.numeric(unique(df$Pos[-grep("M",df$Pos)]))
    names(ids) <- ids
    welch <- sapply(ids, function(i) oneway.test(Av.motif.imp~Base, df[df$Pos==i,], var.equal=FALSE)$p.value) # Welch one-way test
    welch <- p.adjust(welch, method="fdr")
    welch <- ifelse(welch < 0.01, "*", "")
    
    boxplot <- ggplot(df, aes(factor(Pos), Av.motif.imp, fill=Base)) +
      #geom_hline(yintercept = 0, linetype="dashed") +
      geom_boxplot(outlier.size = -1, alpha=0.8, size=0.5) +
      scale_x_discrete(labels=c(-5:-1,
                                strsplit(motif, "")[[1]],
                                paste0("+", 1:5))) +
      ylab(paste0(motif, " importance")) +
      scale_fill_manual(values = col_scheme) +
      theme(panel.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line.y=element_line(colour="black"), axis.line.x=element_blank(),
            axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=14, colour="black"),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            legend.key=element_rect(fill=NA, colour=NA),
            legend.text = element_text(size=11, colour="black"),
            legend.title = element_text(size=12, colour="black"),
            plot.title = element_text(size=16, hjust = 0.5, colour="black"), plot.subtitle = element_text(size=14, hjust = 0.5))
    
    if(class=="dev" & motif=="GATAA"){yl1 <- -0.05; yl2 <- 0.25; yl3 <- 0.26}
    if(class=="dev" & motif=="TGA.TCA"){yl1 <- -0.01; yl2 <- 0.37; yl3 <- 0.375}
    if(class=="dev" & motif=="CATCTG"){yl1 <- -0.02; yl2 <- 0.13; yl3 <- 0.135}
    if(class=="dev" & motif=="GAGAG"){yl1 <- -0.018; yl2 <- 0.078; yl3 <- 0.08}
    if(class=="hk" & motif=="ATCGAT"){yl1 <- -0.007; yl2 <- 0.49; yl3 <- 0.5}
    if(!is.null(yl1)){
      boxplot <- boxplot + geom_text(data=data.frame(x=names(welch), y=yl2), aes(x,y, label=welch), color="black", size=6, inherit.aes = F) +
        coord_cartesian(ylim = c(yl1, yl3))
    }
    yl1 <- NULL
    
    # logo of top quintile
    quintiles <- quantile(df$Av.motif.imp, prob = seq(0, 1, 0.1))
    df <- df[df$Av.motif.imp>=as.numeric(quintiles[length(quintiles)-1]),]
    prop_matrix <- apply(table(df$Base, df$Pos), 2, function(x) x/sum(x))
    # extend 4 positions
    prop_matrix <- cbind(rep(0.25,4), rep(0.25,4), prop_matrix, rep(0.25,4), rep(0.25,4))
    
    plot.logo <- function(x, method='bits'){
      ggseqlogo(x, method=method, seq_type='dna', ncol=1) +
        theme(panel.background = element_rect(fill="white",colour="white"),
              panel.grid = element_blank(), axis.line=element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(),
              axis.title = element_blank(),
              legend.key=element_blank(),
              legend.text = element_text(size=11, colour="black"),
              legend.title = element_text(size=12, colour="black"),
              plot.title = element_text(size=14, hjust = 0.5, colour="black"), plot.subtitle = element_text(size=14, hjust = 0.5),
              legend.position = "right")
    }
    
    print(plot_grid(plot.logo(prop_matrix) + ggtitle(paste0("Logo of 90th percentile of ", motif, " instances (", class, " enh)")),
                    plot_grid(boxplot, heatmap, align = "v", nrow = 2, rel_heights = c(1/4, 3/4)),
                    align = "v", nrow = 2, rel_heights = c(1/7, 6/7)))
    
  }
  dev.off()
  
  print(class)
}


########
# Fig 4D & Fig S16B,D
########

df <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/GATAA_replace_flanks_main_log2FC_table.txt")

# paired Wilcoxon signed rank test
wilcox.test(df$dev_wt[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nstrong by weak" & df$Experiment_detail_flank_length %in% "2bp flanks"],
            df$dev_mut[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nstrong by weak" & df$Experiment_detail_flank_length %in% "2bp flanks"], paired=T)$p.value # **
wilcox.test(df$dev_wt[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nweak by strong" & df$Experiment_detail_flank_length %in% "2bp flanks"],
            df$dev_mut[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nweak by strong" & df$Experiment_detail_flank_length %in% "2bp flanks"], paired=T)$p.value # *

wilcox.test(df$dev_wt[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nstrong by weak" & df$Experiment_detail_flank_length %in% "5bp flanks"],
            df$dev_mut[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nstrong by weak" & df$Experiment_detail_flank_length %in% "5bp flanks"], paired=T)$p.value # ****
wilcox.test(df$dev_wt[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nweak by strong" & df$Experiment_detail_flank_length %in% "5bp flanks"],
            df$dev_mut[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace\nweak by strong" & df$Experiment_detail_flank_length %in% "5bp flanks"], paired=T)$p.value # *

g1_log2FC <- ggplot(df[complete.cases(df$Strongest_log2FC_instance),], aes(Var, dev_log2FC_wt_mut, fill=Var)) +
  geom_boxplot(outlier.size = -1) +
  facet_grid(~Experiment_detail_flank_length) +
  scale_y_continuous("log2 FC to wt", breaks = seq(-10,10,1)) +
  coord_cartesian(ylim = c(-4.7,2.8)) +
  xlab("Motif flanks replaced") +
  geom_hline(yintercept = 0, linetype="dashed", col="grey40") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3, "Dark2")[3],
                               RColorBrewer::brewer.pal(3, "Dark2")[1])) +
  guides(fill=F) +
  theme(axis.text.x = element_text(size=12))

ggsave("Fig4D.pdf", g1_log2FC, width = 6.5, height = 4)

## Fig S16D
df <- read.delim("https://data.starklab.org/almeida/DeepSTARR/Figures_data/GATAA_swap_and_mutation_main_table.txt")
df$Var <- factor(df$Var, levels = c("Mut strong", "Replace weak\n& mut strong", "Mut weak", "Replace strong\n& mut weak"))

# paired Wilcoxon signed rank test
wilcox.test(df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Mut strong" & df$Experiment_detail_flank_length %in% "2bp flanks"],
            df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace weak\n& mut strong" & df$Experiment_detail_flank_length %in% "2bp flanks"], paired=T)$p.value # n.s.
wilcox.test(df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Mut weak" & df$Experiment_detail_flank_length %in% "2bp flanks"],
            df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace strong\n& mut weak" & df$Experiment_detail_flank_length %in% "2bp flanks"], paired=T)$p.value # **

wilcox.test(df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Mut strong" & df$Experiment_detail_flank_length %in% "5bp flanks"],
            df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace weak\n& mut strong" & df$Experiment_detail_flank_length %in% "5bp flanks"], paired=T)$p.value # *
wilcox.test(df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Mut weak" & df$Experiment_detail_flank_length %in% "5bp flanks"],
            df$value[complete.cases(df$Strongest_log2FC_instance) & df$Var %in% "Replace strong\n& mut weak" & df$Experiment_detail_flank_length %in% "5bp flanks"], paired=T)$p.value # ***

gg <- ggplot(df[complete.cases(df$Strongest_log2FC_instance) & complete.cases(df$Var),], aes(Var, value, fill=Var, alpha=Var)) +
  geom_boxplot(outlier.size = -1) +
  facet_grid(~Experiment_detail_flank_length) +
  scale_y_continuous("log2 FC to wt", breaks = seq(-10,10,1)) +
  xlab("Enhancer variants") +
  geom_hline(yintercept = 0, linetype="dashed", col="grey40") +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(3, "Dark2")[3],
                               RColorBrewer::brewer.pal(3, "Dark2")[3],
                               RColorBrewer::brewer.pal(3, "Dark2")[1],
                               RColorBrewer::brewer.pal(3, "Dark2")[1])) +
  scale_alpha_manual(values = c(1,0.5,1,0.5)) +
  guides(fill=F, alpha=F) +
  theme(axis.text.x = element_text(size=10))

ggsave("FigS16D.pdf", gg, width = 8.5, height = 4)


########
# Fig 4E
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
pdf("Fig4E.pdf", height = 2.1, width = 13)

i="chr3L_20269495_20269743_+_wt_dCP"
my.logo(twist_contr_scores$dev[[grep(i, names(twist_contr_scores$dev), fixed = T)]]) + ggtitle(paste0(i, " - dev scores"))

dev.off()
