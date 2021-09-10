#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
options("scipen"=100, "digits"=4)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-i", "--input"),  action = "store", type="character", default=NULL,
              help="$TMP/collapsed_frags.bed", metavar="character"),
  make_option(c("-m", "--MM"),  action = "store", type="numeric", default=1,
              help="$UMI_MM", metavar="character"),
  make_option(c("-c", "--core"),  action = "store", type="numeric", default=1,
              help="number of cores", metavar="character"),
  make_option(c("-o", "--out"),  action = "store", type="character", default="average.txt",
              help="$TMP/reads.filtered.3.bed", metavar="character")
)

opt_parser <- OptionParser(
  usage = "%prog [options]",
  option_list=option_list,
  description = "UMI collapsing"
)
arguments <- parse_args(opt_parser, positional_arguments = TRUE)
opt <- arguments$options

#------------
# LIBRARIES
#------------

suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("stringdist"))

# print options
cat("\nRunning UMI collapsing\n")

opt

#------------
# Prepare data
#------------

test_big <- import.bed(opt$input)
test_big$ID <- paste(test_big@seqnames, test_big@ranges@start, test_big@ranges@start+test_big@ranges@width-1, test_big@strand, sep="_")
lvl <- names(sort(table(test_big$ID)))
test_big$ID <- factor(test_big$ID, levels = rev(lvl))
test_big_sorted <- test_big[order(test_big$ID)]

f3 = function(bar,c=1){
  keep <- vector("numeric",length(bar))
  names(keep) = names(bar)
  while (length(bar)>0) {
    rmv = which(stringdist(bar[1],bar,method="hamming",nthread =c)<=opt$MM)
    keep[names(bar)[1]] = length(rmv)
    bar = bar[-rmv]
  }
  return(keep)
}

#------------
# Run
#------------

# Calculate the number of cores
no_cores <- opt$core
# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl, "stringdist")
clusterExport(cl, "f3")
clusterExport(cl, "opt")

results <- unlist(parLapply(cl, split(test_big_sorted$name,test_big_sorted$ID),
                                        function(x){
                                          names(x) <- 1:length(x)
                                          f3(x)
                                        })
)

# Finish
stopCluster(cl)

test_big_sorted$counts <- as.numeric(results)
# print to see top
test_big_sorted
test_big_sorted <- test_big_sorted[test_big_sorted$counts>0]
test_big_sorted

out <- as.data.frame(test_big_sorted)
out$start <- out$start-1
out$name <- paste(out$name, out$score, sep="_")

write.table(out[,c(1:3,6,9,5)], opt$out, sep="\t", row.names = F, col.names = F, quote = F)

sessionInfo()
