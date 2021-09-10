#!/bin/bash

# Stark Lab in house pipeline. Adapted by Bernardo Almeida 2020

################################################################################
# Requirements
################################################################################

# Programs:
# * Jim Kent's utils
# * bedtools
# * grep-overlap
# * score-overlap
# * median.R
# * hyper.R
# * precomp_pvals.R

#update 2020/03/25 We have only one cluster now, clip, therefore we dont need other configurations
  module load stark.grp/0.0.2
  Rexec="singularity run --app Rscript /groups/stark/software-all/singularity_img/singularity.R.with_pckg.simg "


# Files:
# * /groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes

################################################################################
# Set default values
################################################################################

# params to precompute min. starrs-seq count such that a sign. enrichment can be
# reached given the input count and the overall library sizes
# be cautious when increasing these params - this increases with (nearly) MS*MI
MS=500       # precompute for starr-seq fragment count up to MS
MI=200       # precompute for input fragment count up to MI

WIN=500      # window to evaluate local background
DST=0        # min. distance between peak-summits (=peak width)
EXT=0        # extend fragments to that length
Z=1.67       # z-score for conf-ratio
PVL=0.001    # p-value cutoff for bed file
LPVL=0       # loose p-value cutoff during peak calling
ASSEMBLY=dm3 # assembly/chromosome sizes
STRICT=0     # strict peak calling taking only the lowest abundant strand

MINENR=1     # fold enrichment cutoff for peak calling  (candidates must be > MINENR enriched)

PRUNEP=1     # prune peaks with Daniel's code
memory=4G
#tempDir="/scratch-cbe/users/stark/pipelines/"
tempDir="/tmp/"
################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
                echo >&2 "
$(basename $0) - Call peaks

USAGE: $(basename $0) -e <experiment file> -e <background file> [OPTIONS]
 -e     Input file for experiment (BigBed) [required]
 -b     Input file for background/input (BigBed) [required]
 -o     Prefix for output files (.txt and .bed) [default: prefix for experiment]
 -P     P-value cutoff [default: $PVL]
 -M     Min. fold enrichment cutoff [default: $MINENR]
 -W     Window to evaluate local background [default: $WIN]
 -D     Min. distance between peak summits (=peak width) [default $DST]
        0: use median of all input fragments (e.g. for paired-end data)
        int: use user defined summit-distance/peak-width (e.g. for single-end data; e.g. 600)
 -E     Extend fragments to fixed fragment length E [default: $EXT]
        0: use fragment size itself (e.g. for paired-end data)
        int: use constant fragment length int (e.g. for single-end data)
 -Z     Z-score for correcting enrichments [default: $Z]
 -g     Genome assembly (e.g. dm3, hg19) [default: $ASSEMBLY]
 -L     Loose P-value cutoff during peak-calling (OFF: 0) [default: $LPVL]
 -S     Strict peak calling taking the lowest strand coverage [default: $STRICT]
 -R     Prune peaks (remove peaks that are on flanks of other peaks (0/1) [default: $PRUNEP]
 -m     Memory restrictions for shell sort                                 [default: $memory]

NOTES:
A file with chromosome sizes must bed stored here /groups/stark/genomes/chrom/ASSEMBLY.chrom.sizes.
Uses 4 cores for sorting.
"
        exit 1
fi

################################################################################
# Parse input
################################################################################

while getopts "e:b:o:P:M:W:D:E:Z:g:C:L:S:R:m:" o
do
  case "$o" in
    e) EF="$OPTARG";;
    b) BF="$OPTARG";;
    o) OF="$OPTARG";;
    P) PVL="$OPTARG";;
    M) MINENR="$OPTARG";;
    W) WIN="$OPTARG";;
    D) DST="$OPTARG";;
    E) EXT="$OPTARG";;
    Z) Z="$OPTARG";;
    g) ASSEMBLY="$OPTARG";;
    L) LPVL="$OPTARG";;
    S) STRICT="$OPTARG";;
    R) PRUNEP="$OPTARG";;
    m) memory="$OPTARG";;
   \?) exit 1;;
  esac
done

if [ -z "$EF" -o -z "$BF" ]; then
  echo >&2 "ERROR: -e -b are required!"; exit 1
fi

if [ -z "$EF" -o -z "$BF" ]; then
  echo >&2 "ERROR: input files $EF or $BF not found"; exit 1
fi

if [ $(awk -vP=$PVL -vL=$LPVL 'BEGIN{if(L<P){print 1}else{print 0}}') -eq 1 ]; then
  # LPVL must not be smaller than PVL
  LPVL=$PVL
fi
### update 2019.06.26 -- bc language is not installed in the clip cluster!
#if [ $(echo "if ($MINENR < 1) 1 else 0" | bc) -eq 1 ]; then
#if [ $MINENR < 1 ]; then
#  echo "Warning: Minimal inreachment cant be less than one! so, have changed it to 1"
  # MINENR must be > 1
#  MINENR=1
#fi


################################################################################
# Set chromosome sizes
################################################################################

# Set $INDEX and $SIZES
if [ "$ASSEMBLY" = "dm3" ]; then
  INDEX=/groups/stark/indices/bowtie/dm3_no_chrU_chrUextra_chrM/dm3
  CS=/groups/stark/genomes/chrom/dm3.chrom.sizes
elif [ "$ASSEMBLY" = "ecoli" ]; then
  INDEX=/groups/stark/gerlach/data/ecoli/ecoli
  CS=/groups/stark/gerlach/data/ecoli/ecoli.chrom.sizes
else
  INDEX=/groups/stark/indices/bowtie/${ASSEMBLY}/${ASSEMBLY}
  CS=/groups/stark/genomes/chrom/${ASSEMBLY}.chrom.sizes
fi

# Throw error message if chromosome sizes file does not exist
[ -e "$CS" ] || echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"


################################################################################
# Set output file name, TMP directory, compute number of fragments, p-vals
################################################################################

## Path in which the script is located
DIR="$(cd "$(dirname "$0")" && pwd)"

## If output prefix is not set, use prefix of EF
if [ -z "$OF" ]; then
  OF=${EF%.bb}
else
  if [ ! -d "$(dirname ${OF})" ]; then
    mkdir -p "$(dirname ${OF})"
  fi
fi

## Temporary directory
TMP=$(mktemp -d -p "${tempDir}")
trap "rm -rf $TMP" EXIT

## Get total number of fragments in each file
E=$(bigBedInfo $EF | awk '(/^itemCount/){gsub(/,/,"",$2);print $2}')
B=$(bigBedInfo $BF | awk '(/^itemCount/){gsub(/,/,"",$2);print $2}')

## Pre-compute cutoff values for sign. enrichment given STARR-seq coverage values 0 to MS
## (only try input coverage values 0 to MI)
## Format is STARR, Input, .., ..
 precomp_pvals.R -MS $MS -MI $MI -LS $E -LI $B -p $LPVL -e $MINENR -FU H > $TMP/pvalues

## MINCOV : min. STARR-seq coverage for peak calling (candidates require >= MINCOV fragments)
MINCOV=$(awk '{if($2==1){M=$1}} END{if(M<1){M=1} print M}' $TMP/pvalues)

################################################################################
# Extend reads if necessary
################################################################################

if [ "$EXT" = "0" ]; then
  # convert to bed format, do NOT extend
  bigBedToBed $EF stdout > $TMP/EF.bed
  bigBedToBed $BF stdout > $TMP/BF.bed
  if [ "$DST" = "0" ]; then
  # get median fragment length to use as distance D between summits
    DST=$(awk '{print $3-$2}' $TMP/BF.bed | median.R -i - | awk '{print int(($1/10)+0.5)*10}')
  fi
else
  # do extend fragment size to a total size of EXT
  # use median fragment length (here: EXT) as distance D between summits
  if [ "$DST" = "0" ]; then
    DST=$EXT
  fi
  # redefine EXT to be the extension
  EXT=$(bigBedToBed $EF stdout | head -n 1 | awk -v E=$EXT '{print E-($3-$2)}')
  # convert to bed format, extend
  bigBedToBed $EF stdout | \
      awk -v OFS="\t" -v E=$EXT -v C=$CS 'BEGIN{ while(getline<C){chr[$1]=$2} } { if($6=="+"){$3=$3+E;if($3>chr[$1]){$3=chr[$1]}}else{$2=$2-E;if($2<0){$2=0}} print $0 }' | \
      sort -k1,1 -k2,2n -k6,6 -k5,5n -k4,4 -S $memory --parallel 4 > $TMP/EF.bed
  bigBedToBed $BF stdout | \
      awk -v OFS="\t" -v E=$EXT -v C=$CS 'BEGIN{ while(getline<C){chr[$1]=$2} } { if($6=="+"){$3=$3+E;if($3>chr[$1]){$3=chr[$1]}}else{$2=$2-E;if($2<0){$2=0}} print $0 }' | \
      sort -k1,1 -k2,2n -k6,6 -k5,5n -k4,4 -S $memory --parallel 4 > $TMP/BF.bed
fi



################################################################################
# Get coverage for Input & summit candidates for STARR-seq (from coverage)
################################################################################

if [ "$STRICT" = "1" ]; then
    ## strict counting - STARR-seq counts are corrected UP
    ## to the level of the strand with the least depletion

    ## get strand-specific coverage files for EF and non-strand-specific for BF (bedgraph format)
    PIDS=""
    for s in "+" "-"; do
	cat $TMP/EF.bed | genomeCoverageBed -i stdin -bg -g $CS -strand "$s" > $TMP/EF.coverage.$s &
	PIDS="$PIDS $! "
    done
    cat $TMP/BF.bed | genomeCoverageBed -i stdin -bg -g $CS > $TMP/BF.coverage.b &
    PIDS="$PIDS $! "

    wait $PIDS

    ## merge strand-specific coverage files to one overall file
    ## -> record only region centers as summit candidates (if they pass a coverage cutoff)
    ## (still missing: merge touching regions with identical fragment count)
    bedtools unionbedg -i $TMP/EF.coverage.+ $TMP/EF.coverage.- | \
        awk -v OFS="\t" '{p=int(($2+1+$3)/2); S=($4<=$5 ? $4*2 : $5*2);  print $1,p,p,S}' > $TMP/summit-candidates
    ## clean up strand-specific coverage files
    for s in "+" "-"; do
	rm $TMP/EF.coverage.$s
    done

else
    ## normal counting - sum of both strands is used, we obtain this directly...
    ## for STARR-seq, only region centers as summit candidates are recorded (if they pass a coverage cutoff)

    PIDS=""
    cat $TMP/BF.bed | genomeCoverageBed -i stdin -bg -g $CS > $TMP/BF.coverage.b &
    PIDS="$PIDS $! "
    cat $TMP/EF.bed | genomeCoverageBed -i stdin -bg -g $CS | awk -vM=$MINCOV '($4>=M){p=int(($2+1+$3)/2); print $1,p,p,$4}' > $TMP/summit-candidates &
    PIDS="$PIDS $! "

    wait $PIDS

fi


################################################################################
# Filter summit candidates based on their enrichment over input
################################################################################

## this first intersects the 1 nts long summit candidates with the input coverage as it is much faster
## filtering is based on whether the enrichment and pvalue (precomp.) cutoffs can be reached
## in a 2nd step, summit regions are expandedto +/- DST and a non-overlapping set of candidates is computed using score-overlap

awk -vCS=$CS 'BEGIN{while((getline<CS)>0){C[$1]=$2}} {if($1!=OC){OC=$1; print $1,1,C[$1],0} print $1,$2+1,$3,$4}' $TMP/BF.coverage.b | \
    grep-overlap -s ' ' -c 0 $TMP/summit-candidates - | awk '{k=($1" "$2); if(k!=o){ if(o!=""){if(b<1){b=1} print o,e,b } o=k } e=$4; b=$8 } END{if(b<1){b=1} print o,e,b }' | \
    awk -vPVS=$TMP/pvalues -vME=$MINENR -vE=$E -vB=$B -vD=$DST -vW=$WIN -vCS=$CS 'BEGIN{R=E/B; DE=1/ME; while((getline<PVS)>0){MINSTARR[$2]=$1} MAX=$2; while((getline<CS)>0){C[$1]=$2}}
         { P=0; if($4 in MINSTARR){ if($3>=MINSTARR[$4]){P=1} }else{ if(($3>=MINSTARR[MAX]) && (R*$4/$3<DE)){P=1} }
           if(P){ s=$2-int(D/2); if(s<1){s=1} b=$2-W; if(b<1){b=1} e=$2+W; if(e>C[$1]){e=C[$1]} print $3,$1,s,$2+int(D/2),$1,b,e,$2,$3 } }' | \
    score-overlap > $TMP/summits.b


################################################################################
# Get fragment counts for summit positions
################################################################################

# Below, I use as input the max(local input at position OR average of a W window)
# We don't use a pseudocount but set input to 1 if it is 0
# Future: use rpy to calculate binom. p-value inside python -> should be faster

cat $TMP/BF.coverage.b | \
  awk -vCS=$CS 'BEGIN{while((getline<CS)>0){C[$1]=$2}} {if($1!=OC){OC=$1; print $1,1,C[$1],0} print $1,$2+1,$3,$4}' | \
  grep-overlap -s ' ' -c 0 $TMP/summits.b - | \
  awk '{k=($1" "$4);if(k!=o){if(o!=""){print o,e,b,int(s/l+0.5);s=0}o=k;l=$3-$2+1} if($9>0){s+=$12*$9} if($4>=$7 && $4<=$8){e=$5;b=$9}} END{print o,e,b,int(s/l+0.5)}' | \
  awk -vE=$E -vB=$B -vME=$MINENR 'BEGIN{R=E/B} {b=$4>$5?$4:$5; if(b<1){b=1} if(($3/b)/R>ME){print $3,E,b,B,$3,$3+b,E,E+B,$1,$2,$3,b,E,B,($3/b)/R}}' | \
  ${DIR}/conf_ratio.py -z $Z | hyper.R -i - -m 5 -n 6 -M 7 -N 8 -o TRUE -u FALSE | \
  cut -f 9-1000 | awk -vP=$LPVL '$9<=P' | \
  sort -k1,1 -k2,2n -S $memory --parallel 4 > $TMP/peaks.b


################################################################################
# Peak pruning
# Remove peaks which are on the flanks of a peak primary peak.
# Highest fragment count does not overlap with the summit position but is
# shifted towards the ends.
################################################################################

# Print all full peaks regions, intersect with fragments, count height for each position
# in the peak, report diference of original summit and new highest point,
# sort by highest point and then by difference to summit position,
# unique peak chr and old summit, and report the difference in bp, in %, and the new height
# filter on more than 10% different than original summit

if [ "$PRUNEP" = "1" ]; then
  cat $TMP/peaks.b | \
    awk -vDST=$DST '{D=DST/2} {for (i=-D;i<=D;i++) {print $1,$2+i,$2+i,$2,$0}}' | \
    grep-overlap -c 0 -s " " - <(awk '{$2+=1; print $0}' $TMP/EF.bed) | \
    awk '{x[$1" "$2" "$3" "$4" "$7" "$8" "$9" "$10" "$11" "$12" "$13]++} END{for (i in x) {print x[i],i}}' | \
    awk 'function abs(x){return (((x < 0) ? -x : x) + 0)} {print abs($5-$3),$0}' | \
    sort -k2,2nr -k1,1n -S $memory --parallel 4 | \
    awk -vD=$DST -vOFS="\t" '!x[$3" "$6]++{print $3,$6,$7,$8,$9,$10,$11,$12,$13,$1,$1/D*100,$2}' | \
    sort -k8,8gr -S $memory --parallel 4 | \
    awk '$11<=10' | \
    cut -f1-9 > $OF.peaks.txt
else
  sort -k8,8gr -S $memory --parallel 4 $TMP/peaks.b > $OF.peaks.txt
fi


################################################################################
# Filter peaks and write a BED file
################################################################################

cat $OF.peaks.txt | \
  awk -vD=$DST -vCS=$CS -vP=$PVL -vOFS="\t" 'BEGIN{while((getline<CS)>0){S[$1]=$2}}
      ($9<=P){s=$2-int(D/2)-1;if(s<0)s=0;e=$2+int(D/2);if(e>S[$1])e=S[$1];print $1,s,e,"peak_"NR,$3<=1000?$3:1000,".",($2-1>=0?$2-1:0),$2,"50,0,0"}' | \
  sort -k1,1 -k2,2n -S $memory --parallel 4 > $OF.peaks.bed

################################################################################
# Convert BED to BigBed
################################################################################

bedToBigBed $OF.peaks.bed $CS $OF.peaks.bb 2>/dev/null

# Exit
rm -rf ${TMP}

exit 0
