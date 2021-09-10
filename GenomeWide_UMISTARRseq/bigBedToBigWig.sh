#!/bin/bash
set -o errexit
set -o pipefail

# Stark Lab in house pipeline. Adapted by Bernardo Almeida 2020

################################################################################
# Requirements
################################################################################

# Programs:
# * Jim Kent's utils (http://hgwdev.cse.ucsc.edu/~kent/src/)
# * bedtools (http://code.google.com/p/bedtools/)

module load stark.grp/0.0.2

################################################################################
# Set default values
################################################################################

ASSEMBLY="dm3"  # Genome assembly
OUTFILE=""      # Outfile, if not supplied ${input%.bb}.bw
EXT="0"         # Fragment length of the experiment
EXT2="0"        # Addititional extension
NORM="1"        # Normalize read coverge to 1 million mapped reads
TYPE="libSize" # Normalization type
SCALING_FILE="" # Normalize to the precalculated scalling factor
STRAND="0" #make strand-specific tracks
add="none"

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - Compute a BigWig file from a BigBed file of mapped reads

USAGE: $(basename $0) -i <BigBed input file> [OPTIONS]
 -i     Input file (BigBed)                                                                 [ required ]
 -o     Output file (BigWig)                                                                [default: <input file>.bw ]
 -g     Genome assembly (e.g. dm3, hg19)                                                    [default: $ASSEMBLY]
 -e     Extend mapped read coordinates to total length of n (off: 0)                        [default: $EXT ]
        Set to 0 for paired-end data
        Set to fragment size (e.g. 600) for single-end data to extend fragment coordinates
 -E     Additional extension                                                                [default: $EXT2 ]
        Enables the artificial extension of reads.
 -t     Type of normalization: library size [libSize], pre-calculated scaling factor[scaling], or
        none                                                                                [default: $TYPE ]
 -n     Normalize the coverage to 1 million mapped reads (0/1)                              [default: $NORM]
 -F     File with precalculated scalling factor.                                            [default: none ]
        File should be tabular separated and contain two column: experimentId and its
        scaling factor. Experiment ID should be exactly the same as in experiment file,
        usually it is a field Outfile
 -s     Make strand-specific bigWig files [0/1]                                             [default: 0 ]
 -a     'Add' to the fragment start. Used for PRO-seq, where we need only start postion of  [default: $add]
        the fragment ( -a 1)

NOTES:
The program requires a file with the chromosome sizes of your assembly stored
here /groups/stark/genomes/chrom/ASSEMBLY.chrom.sizes.

The program uses up to 4 CPU cores, and about 2GB of RAM for a 100MB BigBed file.

Also script required
"
  exit 1
fi

################################################################################
# Parse input
################################################################################

while getopts "i:o:g:e:n:E:F:t:s:a:" o
do
  case "$o" in
    i) INFILE="$OPTARG";;
    o) OUTFILE="$OPTARG";;
    g) ASSEMBLY="$OPTARG";;
    e) EXT="$OPTARG";;
    E) EXT2="$OPTARG";;
    t) TYPE="$OPTARG";;
    F) SCALING_FILE="$OPTARG";;
    n) NORM="$OPTARG";;
    s) STRAND="$OPTARG";;
    a) add="$OPTARG";;
   \?) exit 1;;
  esac
done

################################################################################
# Set chromosome size file
################################################################################

# Set $INDEX and $SIZES
if [ "$ASSEMBLY" = "dm3" ]; then
  INDEX=dm3_no_chrU_chrUextra_chrM/dm3
  SIZES=dm3.chrom.sizes
else
  INDEX=${ASSEMBLY}/${ASSEMBLY}
  SIZES=${ASSEMBLY}.chrom.sizes
fi

# Throw error message if chromosome sizes file does not exist
if [ ! -e "$SIZES" ]; then
  echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"
  exit 1
fi


#some patch to be able to run this script for old data sets with old parameters
#TYPE parameter has priority on the NORM
if [ "$TYPE"=="none" ]; then
  NORM="0"
fi

#but if it is not set, then NORM will set priority
if [ -z "$TYPE" ]; then
 if [ "$NORM" = "1" ]; then
   TYPE="libSize"
 else
   TYPE="none"
 fi
fi



################################################################################
# Run main program
################################################################################

# Get number of total fragments

N=$(bigBedInfo $INFILE | awk '(/^itemCount/){gsub(/,/,"",$2);print $2}')

#Get experiment ID from the file name

expId=$(echo $INFILE | perl -ne 'chomp($_); ($id)=$_=~/\/*([^\/.]+)\.[^\/]*bb/; print $id;')

if [ $expId == "" ]; then
 echo "Experiment is not found!\n"
 exit 1
fi

#get scalling factor
if [ "$TYPE" = "scaling" ]; then
  if [ -z "$SCALING_FILE" ]; then
   echo "Scaling file is missing! "
   exit 1;
  fi
  
  SFline=$(grep -w $expId $SCALING_FILE)
  if [ "$SFline" == "" ]; then
   echo "cant find scalling factor for sample $expId!"
   exit 1;
  fi

  SF=$(grep -w $expId $SCALING_FILE| awk '{print $2}')

fi

echo "$SF"
echo "normalization..."

# Get mapped read length and compute extension paramter E
# Set to zero if EXT should be zero
E=$(bigBedToBed $INFILE stdout -maxItems=1 | awk -v E=$EXT 'NR==1{if(E==0){print "0"}else{print E-($3-$2)}}')

# Error if extension lower than 0
if [ "$E" -lt "0" ]; then
  echo >&2 "ERROR: Total fragment length should be larger than read length!"
  exit 1
fi

# Set outfile to infile.bw or to the user given name
if [ "$OUTFILE" = "" ]; then
  OUTFILE=${INFILE%.bb}.bw
else
  OUTFILE=$OUTFILE
fi


# Compute coverage and get BedGraph file
# Convert BedGraph file to BigBed
# If required extend the reads prior to computing the coverage
# EXT2 is and additional extension parameter

#update 21.01.19 - sometimes we need strand-specific tracks

if [ "$STRAND" = "1" ]; then
  strandList=("+" "-")
else
  strandList="none"
fi

for strand in  ${strandList[*]}; do
  echo $strand

   if [ "$strand" = "+" ]; then
    OUTFILE=${OUTFILE%.bw}_ps.bw
    if [ $add != "none" ]; then
       bigBedToBed $INFILE stdout |  awk -vOFS="\t" -vadd=$add '($6=="+"){$3=$2+add; print}' > ${OUTFILE}.bed
    else
      bigBedToBed $INFILE stdout |  awk -vOFS="\t"  '($6=="+"){print}' > ${OUTFILE}.bed
    fi
  elif [ "$strand" = "-" ]; then
    OUTFILE=${OUTFILE%_ps.bw}_ns.bw
      if [ $add != "none" ]; then
       bigBedToBed $INFILE stdout |  awk -vOFS="\t" -vadd=$add '($6=="-"){$2=$3-add; print}' > ${OUTFILE}.bed
     else
        bigBedToBed $INFILE stdout |  awk -vOFS="\t"  '($6=="-"){ print}' > ${OUTFILE}.bed
     fi
 else
    bigBedToBed $INFILE ${OUTFILE}.bed
 fi

  echo $OUTFILE
cat ${OUTFILE}.bed |\
  if [ "$E" = "0" -a "$EXT2" = "0" ]; then
    cat
  else
    awk -vE=$E -vE2=$EXT2 -vC=$SIZES -vOFS="\t" '
      BEGIN {while(getline<C){chr[$1]=$2}}
      {
        if($6=="+"){$3=$3+E+E2}else{$2=$2-E-E2}
        if($2<0){$2=0}
        if($3>chr[$1]){$3=chr[$1]}
        print $0
      }'
  fi |\
  genomeCoverageBed -i stdin -bg -g $SIZES |\
  if [ "$TYPE" = "libSize" ]; then
       awk -vN=$N -vOFS="\t" '{print $1,$2,$3,1e6*$4/N}'
  elif [ "$TYPE" = "scaling" ]; then    #scaling normalization
    awk -vsf=$SF -vOFS="\t" '{print $1,$2,$3,$4*sf}'
  else   #no normalization
    awk -vN=$N -vOFS="\t" '{print $1,$2,$3,$4}'
  fi | \
  if [ "$strand" = "-" ]; then
    awk -vOFS="\t"  '{print $1,$2,$3,-1*$4}'
  else
    awk -vOFS="\t"  '{print $1,$2,$3,$4}'
  fi | \

  wigToBigWig stdin $SIZES $OUTFILE
  rm ${OUTFILE}.bed
done


# Exit
exit 0
