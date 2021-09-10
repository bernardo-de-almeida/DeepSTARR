#!/bin/bash
set -o errexit
set -o pipefail

# Stark Lab in house pipeline. Adapted by Bernardo Almeida 2020

################################################################################
# Requirements
################################################################################

# Programs:
# * samtools
# * bowtie
# * slippage_filter_pe.sh
# * Jim Kent's utils
# * bedtools

OS_id=$(echo $LMOD_SYSHOST)

if [ "$OS_id" == "" ]; then
  module load samtools/0.1.18
  module load kent-ucsc/3.8f6f5e0a1cb75
  module load bowtie/0.12.9
  module load bedtools/2.19.1
elif [ "$OS_id" == "IMPIMBA-2" ]; then  #ii2 environment
  module load samtools/0.1.19-foss-2017a
  module load ucsc-kent-utils/356-foss-2017a
  module load bowtie/1.2.2_p1-foss-2017a
  module load bedtools/2.27.1-foss-2017a
  module load gcc/6.3.0-2.27
  Rexec="Rscript "
else  #new clip environment
  module load stark.grp/0.0.2
  module load cutadapt/1.18-foss-2018b-python-3.6.6
  Rexec="Rscript  "
fi

# Notes:
# Needs at least 4 cores

################################################################################
# Set default values
################################################################################

ASSEMBLY="dm3"
BARCODES14="*" #this barcode is used for UMI sequences as well, so by default it is better to have *
BARCODES12="*"
BARCODE14_LEN=6 # Trueseq IDX ("reverse") BC
BARCODE12_LEN=8 # Nextera i5 ("forward") BC
CUT="0"
CUT2="1000"
FILTER="A"
KEEPSEQ="0"
KEEPMM="0"
KEEPLOG="1"
GAP="2000"
BOWTIE_MM="3"
FASTA="0"
USE=4
KEEP="0"
UMI_flag="FALSE"
UMI_MM=1
OUTDIR="data/"
VERBOSE="FALSE"
tempDir="/tmp/"

################################################################################
# Help
################################################################################

function usage {
echo "
$(basename $0) - Extract barcoded reads from a paired-end BAM file and map the reads against a genome assembly

USAGE: $(basename $0) -i <BAM input file> -o <BigBed output file> -b <barcodes> [OPTIONS]

## INPUT/OUTPUT ##
 -i|--input      Input file (BAM) or gzipped fasta                                            [required]
 -o|--output     Output file without path (output format: BigBed)                             [required]
 --outdir        Output folder (output format: BigBed)                                        [default: data/]


 -F              Save gzipped FASTA and exit without mapping (0/1)
                 Reads will be saved as outfile_1.fasta.gz and outfile_2.fasta.gz             [default: $FASTA]
 -f              Filters                                                                      [default: $FILTER]
                      (0) Report all mapped reads
                      (1) Report fragments collapsed on identical chromsome, start, end,
                              and strand information
                      (2) Report unique fragments passing the sequence-slippage filter
                              for paired-end read clusters with same start or same end
                      (3) Report fragments collapsed on single start per chromsome, end, strand,
                              and single end per chromosome, start, strand information
                      (A) Report both all mapped reads AND unique fragments after slippage-filter,
                              i.e. report results from both filters 0 and 2

 -k              Keep intermediate files, for debugging purpuses (0/1)                        [default: $KEEP]


 -v              Verbose, to print debugging information                                      [default: none ]


 ## DATA options ##

 -b             Barcodes applied to column 14                                                     [required]
                (separated by |, enclosed in \" \"; use * as wildcard)
                Note: you can not pass * as a command line option.
                If you have an old bam file with the barcode in column 12 and not in 14,
                you must leave barcodes 14 empty and provide barcode 12. Barcodes 14 will be set to "*" internally.
 -B             Barcodes applied to column 12                                                     [default: $BARCODES12]
                (separated by |, enclosed in \" \"; use * as wildcard)
 -l             Barcode length for BC in column 14                                                [default: $BARCODE14_LEN]
 -L             Barcode length for BC in column 12                                                [default: $BARCODE12_LEN]
 --umi          Experiment was prepared using UMI labeling                                        [default: $UMI_flag]
                The reads will be collapsed based on position and UMI sequence identity
 -g             Genome assembly (e.g. dm3, hg19)                                                      [default: $ASSEMBLY]

 ## MAPPING options ##

 -G             Maximum gap size between paired-end reads                                             [default: $GAP]
 -c             Cut reads that are mapped to n base pairs (off: 0)                                    [default: $CUT]
 -C             Do not report mapped reads, but reads cut to n base pairs                             [default: $CUT2]
                (applies only if -c > 0; n > readlength won't cut)
 -s             Report read/fragment sequences tags (0/1)                                             [default: $KEEPSEQ]
 -m             Maximum number of mismatches per indidual read (0/1/2/3)                              [default: $BOWTIE_MM]
 --UMI_mm       Maximum number of mismatches in UMI sequnece (0/1/2/3)                                [default: $UMI_MM]

 -M             Report number of mismatches (column 5 of BED/BB file; 0/1)                            [default: $KEEPMM]
 -u             Number of processes to use (will be overwritten by NSLOTS on SGE, see below)          [default: $USE]
 -r             Report mapping statistics (0/1)                                                       [default: $KEEPLOG]


NOTES:
The program will use \$NSLOTS cores (set by SGE) or by user-defined number (-u option; minimum of 4 cores required).
In addition the programs bowtie and slippage_filter_pe.sh
must be installed on the machine. Note that you need the bowtie index files for your genome
assembly.

The program creates about 10GB of temporary files in /tmp (35 million 50mer reads, paired-end
data).
"
exit 0
}


################################################################################
# Parse input and check for errors
################################################################################

if [[ ! $@ ]]; then
    usage
fi

ARGS=`getopt -o "i:o:b:B:l:L:g:G:c:C:s:m:f:F:u:r:k:hp:vp" -l "input:,exp:,output:,i5:,idx:,i5_len:,idx_len:,genome:,gap:,cut:,cut2:,filter:,keepseq:,keeplog:,cores:,umi,outdir:,help" \
      -n "$0" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
  then
    exit 1
  fi

# set parameters to preprocessed string $ARGS
eval set -- "$ARGS"

#echo "$ARGS"

while true;
do
  case "$1" in
    -i|--input)
      if [ -n "$2" ]; then
        INFILE=$2
      fi
      shift 2 ;;
    -o|--output)
     if [ -n "$2" ]; then
      OUTFILE=$2
     fi
     shift 2 ;;
    --outdir)
      if [ -n "$2" ]; then
       OUTDIR=$2
      fi
      shift 2 ;;
   -b|--idx)
      if [ -n "$2" ]; then
       BARCODES14=$2
      fi
      shift 2 ;;
   -B|--i5)
      if [ -n "$2" ]; then
       BARCODES12=$2
      fi
      shift 2 ;;
   -l|--idx_len)
      if [ -n "$2" ]; then
       BARCODE14_LEN=$2
      fi
      shift 2 ;;
   -L|--i5_len)
       if [ -n "$2" ]; then
        BARCODE12_LEN=$2
       fi
       shift 2 ;;
    -g|--genome)
        if [ -n "$2" ]; then
           ASSEMBLY=$2
        fi
        shift 2 ;;
    -G|--gap)
        if [ -n "$2" ]; then
            GAP=$2
        fi
        shift 2 ;;
    -c|--cut)
        if [ -n "$2" ]; then
          CUT=$2
        fi
        shift 2 ;;
    -C|--cut2)
        if [ -n "$2" ]; then
            CUT2=$2
        fi
        shift 2 ;;
    -s|--keepseq)
        if [ -n "$2" ]; then
         KEEPSEQ=$2
        fi
        shift 2 ;;
    -m)
       if [ -n "$2" ]; then
        BOWTIE_MM=$2
       fi
       shift 2 ;;
   -f|--filter)
       if [ -n "$2" ]; then
           FILTER=$2
       fi
       shift 2 ;;
    -F)
       if [ -n "$2" ]; then
        FASTA=$2
       fi
      shift 2 ;;
    -u|--cores)
        if [ -n "$2" ]; then
          USE=$2
        fi
        shift 2 ;;
    -r|--keeplog)
        if [ -n "$2" ]; then
          KEEPLOG=$2
        fi
        shift 2 ;;

    -k)
      if [ -n "$2" ]; then
        KEEP=$2
      fi
      shift 2 ;;

     --umi)
           UMI_flag="TRUE"
          shift ;;

      -v)
       VERBOSE="TRUE"
       shift ;;
       -h|--help)
           usage
           shift ;;

      --)
           shift
           break ;;
    esac
done


# while getopts "i:o:b:B:l:L:g:G:c:C:s:m:f:F:u:r:" o
# do
#   case "$o" in
#     i) INFILE="$OPTARG";;
#     o) OUTFILE="$OPTARG";;
#     b) BARCODES14="$OPTARG";;
#     B) BARCODES12="$OPTARG";;
#     l) BARCODE14_LEN="$OPTARG";;
#     L) BARCODE12_LEN="$OPTARG";;
#     g) ASSEMBLY="$OPTARG";;
#     G) GAP="$OPTARG";;
#     c) CUT="$OPTARG";;
#     C) CUT2="$OPTARG";;
#     s) KEEPSEQ="$OPTARG";;
#     m) BOWTIE_MM="$OPTARG";;
#     f) FILTER="$OPTARG";;
#     F) FASTA="$OPTARG";;
#     u) USE="$OPTARG";;
#     r) KEEPLOG="$OPTARG";;
#    \?) exit 1;;
#   esac
# done

# note you can not pass * as a command line option. If you have an old bam file with the
# barcode in column 12 you must leave barcodes 14 empty and provide barcode 12.
# Barcodes 14 will be set to * if it is empty and Barcodes 12 greater then 1

##############################################
# check mandatory input and print variables stdout
##############################################


#if [ ${BARCODES14} -eq 0 -a ${BARCODES12} -gt 1 ];then
#  BARCODES14="*"
#fi

if [ ${BARCODES14} == "UMI" ]; then
  BARCODES14="*"
fi

if [ -z "$INFILE" -o -z "$OUTFILE" -o -z "$BARCODES14" ]; then
  echo >&2 "ERROR: -i -o -b are required!"; exit 1
fi

###################
# PRINT USER specified parameters
#TODO - print all parameters to log file
echo $VERBOSE

if [ $VERBOSE = "TRUE" ]; then
eval "for i in {1..30};do printf \"#\";done"
printf "\n"
echo "Script bowtie_pe will run with following parametrs"
eval "for i in {1..30};do printf \"#\";done"
echo " "

echo "Input bam file $INFILE"
echo "Output file prefix $OUTFILE"
echo "Output folder $OUTDIR"
echo "Barcode 14, idx '$BARCODES14', length $BARCODE14_LEN"
echo "Barcode 12, i5 '$BARCODES12', length $BARCODE12_LEN"
if [ $UMI_flag = "TRUE" ]; then
 echo "Reads will be collapsed by UMI, number of mismatches $UMI_MM "
fi

echo " "
echo "Genome assembly to map data $ASSEMBLY"
echo "Number of mismatches, bowtie mapping: $BOWTIE_MM"

eval "for i in {1..30};do printf \"#\";done"
echo " "
echo "Additional parameters"
eval "for i in {1..30};do printf \"#\";done"
echo " "
echo "Tmp folder to process data ${tempDir}"
echo "Keep intermediate files: $KEEP"
echo "Get fastq files only: $FASTQ"
echo " "
fi


###### flags

if [ $UMI_flag = "TRUE" ]; then
  KEEPSEQ=0
fi

################################################################################
# Set index and chromosome sizes
################################################################################

# Set $INDEX and $SIZES
if [ "$FASTA" = "0" ]; then

  if [ "$ASSEMBLY" = "dm3" ]; then
    INDEX=dm3_no_chrU_chrUextra_chrM/dm3
    SIZES=dm3.chrom.sizes
  else
    INDEX=${ASSEMBLY}/${ASSEMBLY}
    SIZES=${ASSEMBLY}.chrom.sizes
  fi

  # Throw error message if either index or chromosome sizes file does not exist
  if [ ! -e "$SIZES" ]; then
      echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"
      exit -1
  fi
  if [ ! -e "${INDEX}.1.ebwt" ]; then
      echo >&2 "ERROR: No bowtie index files found for genome assembly ${ASSEMBLY}!"
      exit -1
  fi
fi

################################################################################
# Set processor usage
# Either assigned by SGE or USE
# (at least 4 as we anyway need 4 later on)
################################################################################

if [ ! -z "$NSLOTS" ]; then
  USE=$NSLOTS
fi


################################################################################
# Create tmp directory and start log file
################################################################################

# A temporary directory
TMP=$(mktemp -d -p "${tempDir}")
trap "rm -rf $TMP" EXIT

# logs
(
  echo -n "Start: "
  date
  echo "bowtie_pe.sh version XXX"
  echo "executing slippage filter"
  echo "$0 $@"
  echo
) > $TMP/LOG.log

echo $TMP

################################################################################################
# check output file name and if it contain folder name which does not exists, produce an error
################################################################################################
folder=$( echo $OUTFILE | perl -ne 'use File::Basename; chomp($_); $dir=dirname($_); if(! defined $dir){$dir="NA"; } print $dir."\n";' )
fileName=$( basename $OUTFILE )

if [ $folder = "NA" ]; then
  OUTFILE=$OUTDIR/$OUTFILE
fi

if [ ! -d  $folder ]; then
     mkdir -p $folder
fi

if [ ! -d $OUTDIR ]; then
  mkdir -p $OUTDIR
fi


####preprocessing step -
### 1. In case if the input is bam file then execute demultiplexing step
### 2. In case if the input is fasta file, then go to the mapping step directly

################################################################################
# Extract barcoded reads
# Reads are saved in $TMP/reads_1.fa and $TMP/reads_2.fa
################################################################################

#demultiplexing and cutting if needed
#start=`date +%s`


if echo $INFILE | grep -i -q "|"; then

  F_READ=$(echo $INFILE | awk -v FS="|" '{print $1}')
  S_READ=$(echo $INFILE | awk -v FS="|" '{print $2}')
  zcat $F_READ > $TMP/reads_1.fa
  zcat $S_READ > $TMP/reads_2.fa
else

echo "extracting reads.."
(echo -n "Extracting read pairs according to barcodes: "; date) >> $TMP/LOG.log

samtools view -@ $USE $INFILE | \
    awk -vB14=$BARCODES14 -vBCL14=$BARCODE14_LEN -vB12=$BARCODES12 -vBCL12=$BARCODE12_LEN -vO=$TMP/reads -vCUT=$CUT -vCUT2=$CUT2 -vUMI_flag=$UMI_flag '
      BEGIN{
            n=split(B14,T,"|");for(i=1;i<=n;i++){BCS14[substr(T[i],1,BCL14)]=1}
            n=split(B12,T,"|");for(i=1;i<=n;i++){BCS12[substr(T[i],1,BCL12)]=1}
      }
      {

       if( ( ("*" in BCS14) || (substr($14,6,BCL14) in BCS14) ) && ( ("*" in BCS12) || (substr($12,6,BCL12) in BCS12) ) ){
            FID=$1;FR=$10
            if(UMI_flag == "TRUE") {BC14=substr($14,6,BCL14)}
            getline
            if($1==FID){
              if(UMI_flag == "TRUE"){ #add UMI sequence to the read name
                if(CUT>0){
                  print ">"substr(FR,1,CUT2) "_" BC14 "\n" substr(FR,1,CUT)  > (O"_1.fa")
                  print ">"substr($10,1,CUT2)"_" BC14 "\n" substr($10,1,CUT) > (O"_2.fa")
                }else{
                    print ">"FR"_"BC14"\n"FR > (O"_1.fa");
                    print ">"$10"_"BC14"\n"$10 > (O"_2.fa")
                }
              }
            else{ #without UMI_flag
                if(CUT>0){
                    print ">" substr(FR,1,CUT2)  "\n" substr(FR,1,CUT)  > (O"_1.fa")
                    print ">" substr($10,1,CUT2) "\n" substr($10,1,CUT) > (O"_2.fa")
                }else{
                    print ">"FR"\n"FR > (O"_1.fa"); print ">"$10"\n"$10 > (O"_2.fa")
                }
              }
            }else{
                print "ERROR: "FID" does not match "$1; exit(1)
            }
       }else{
            getline
       }
      }'
fi


## Check if fasta files exist and contain data
if [ ! -s $TMP/reads_1.fa -o ! -s $TMP/reads_2.fa ]; then
  echo >&2 "temporary files $TMP/reads_1.fa and/or $TMP/reads_2.fa are empty"
  echo >&2 "input files did not contain any reads or matching barcodes"
  exit -1
fi


## Save FASTA file to outpath and exit
if [ "$FASTA" = "1" ]; then
  gzip -c $TMP/reads_1.fa > ${OUTFILE%.bb}_1.fasta.gz
  gzip -c $TMP/reads_2.fa > ${OUTFILE%.bb}_2.fasta.gz
  rm -rf $TMP
  exit 0
fi
#end=`date +%s`
#runtime=$((end-start))

#echo "done in $runtime"



################################################################################
# Run bowtie
# Output is saved in $TMP/reads.bowtie
################################################################################
#start=`date +%s`


(echo -n "Running bowtie: "; date) >> $TMP/LOG.log

# Bowtie
# with option --norc: bowtie will not attempt to align against the reverse-complement reference strand
bowtie -p $USE -f -X $GAP -v $BOWTIE_MM --norc -m 1 --best --strata $INDEX -1 $TMP/reads_1.fa -2 $TMP/reads_2.fa > $TMP/reads.bowtie 2>> $TMP/LOG.log
echo >> $TMP/LOG.log

# delete fasta files
rm $TMP/reads_1.fa $TMP/reads_2.fa

#end=`date +%s`
#runtime=$((end-start))

#echo "bowtie done in $runtime"

################################################################################
# Parse raw bowtie output to BED format and sort
# Sort keys are chr start end strand mismatches sequence
# Output is saved in $TMP/reads.bowtie.bed --  it is fragments between start and end
################################################################################

#start=`date +%s`

(echo -n "Parsing & sorting output: "; date) >> $TMP/LOG.log

awk -vFS="\t" -vOFS="\t" '{if(NR%2==1){
   split($1,T,"/");
   R[T[2]]=T[1];S=((T[2]==1)?"+":"-");
   B=$4;
   if($8==""){FM=0;FQ="-"}
   else{FQ=$8;TS=$8;gsub(/[^,]/,"",TS);
   FM=1+length(TS)}}
                else{
    split($1,T,"/");
    R[T[2]]=T[1];
    if($8==""){BM=0;$8="-"}
    else{TS=$8;gsub(/[^,]/,"",TS);
    BM=1+length(TS)}
    print $3,B,$4+length($5),R[1]"_"R[2],FM+BM,S}
                          }' $TMP/reads.bowtie | \
    sort -k1,1 -k2,2n -k3,3nr -k6,6 -k5,5n -k4,4 -S 1G --parallel $(echo $USE | awk '{print $1-1}') > $TMP/reads.bowtie.bed

# delete raw bowtie output
# rm $TMP/reads.bowtie

#end=`date +%s`
#runtime=$((end-start))

#echo "Bowtie raw output was convirted to bed file in $runtime"

################################################################################
# Apply filters to collapse reads and/or remove information:
#
#   (0) Report all mapped reads
#   (1) Report fragments collapsed on identical chromsome, start, end,
#       and strand information
#   (2) Report unique fragments passing the sequence-slippage filter
#       for paired-end read clusters with same start or same end
#   (3) Report fragments collapsed on single start per chromsome, end, strand,
#       and single end per chromosome, start, strand information
#   (A) Report both all mapped reads AND unique fragments after slippage-filter,
#       i.e. report results from both filters 0 and 2
#
# Output is saved in $TMP/reads.filtered.bed
################################################################################

(echo -n "Applying filters: "; date) >> $TMP/LOG.log

if [ $FILTER = "0" ]; then
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' $TMP/reads.bowtie.bed > $TMP/reads.filtered.bed
elif [ $FILTER = "1" ]; then
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '(!x[$1" "$2" "$3" "$6]++){print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' $TMP/reads.bowtie.bed > $TMP/reads.filtered.bed
elif [ $FILTER = "2" ]; then
    awk '!x[$4]++' $TMP/reads.bowtie.bed | \
        slippage_filter_pe.sh -m 2 -l 10 -o 1 -d 2000 -- | \
        awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' > $TMP/reads.filtered.bed
elif [ $FILTER = "3" ]; then
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{s=($1":"$2":"$6);e=($1":"$3":"$6);if((!(s in S)) && (!(e in E))){
                                                S[s]=1;E[e]=1; print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}
                                             }' $TMP/reads.bowtie.bed > $TMP/reads.filtered.bed
elif [ $FILTER = "A" ]; then
    # run filter 0 -- fragments
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' $TMP/reads.bowtie.bed > $TMP/reads.filtered.0.bed

    #if UMI - collapse by UMI. For slippage filter we will need to remove variable UMI part
    if [ $UMI_flag = "TRUE" ]; then
      #collapse by UMI and position
      awk -vOFS="\t" '{split($4,T,"_"); print $1,$2,$3,T[2],$5,$6}' $TMP/reads.bowtie.bed|\
      sort -k1,1 -k2,2n -k3,3nr -k6,6 -k4,4 | \
      uniq -c | awk -vOFS='\t' '{print $2,$3,$4,$5,$1,$7}' | \
      sort -k1,1 -k2,2n -k3,3nr -k6,6 -k5,5nr > $TMP/collapsed_frags.bed

      #UMIs with up to UMI_MM mismatches should be considered as same read
      $Rexec STARRseq_UMI_collapsing.R -i $TMP/collapsed_frags.bed -m $UMI_MM -c $USE -o $TMP/reads.filtered.3_unsorted.bed
      #save output from UMI collapsing R script
      #cp $TMP/reads.filtered.3_unsorted.bed  $OUTFILE.reads.UMI.Rscript.bed
      #sort bed file because then I will need to covert it to bigbed
      sort -k1,1 -k2,2n $TMP/reads.filtered.3_unsorted.bed > $TMP/reads.filtered.3.bed

      #slippage filter -- fragments

      awk -vOFS="\t" '{split($4,T,"_"); print $1,$2,$3,T[1]"_"T[3],$5,$6}' $TMP/reads.bowtie.bed | awk '!x[$4]++' | \
          slippage_filter_pe.sh  -m 2 -l 10 -o 1 -d 2000 -- | \
          awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' > $TMP/reads.filtered.2.bed

    else

    # run filter 2
     awk '!x[$4]++' $TMP/reads.bowtie.bed | \
     slippage_filter_pe.sh -m 2 -l 10 -o 1 -d 2000 -- | \
     awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' > $TMP/reads.filtered.2.bed


    fi
fi

# Clean up
rm $TMP/reads.bowtie.bed


################################################################################
# Convert output bed to bigbed format and copy output
################################################################################
echo "convert to bigbed format.."
start=`date +%s`


if [ $FILTER = "A" ]; then

    (echo -n "Running bedToBigBed for filter A-0,reporting all reads: "; date) >> $TMP/LOG.log
    bedToBigBed $TMP/reads.filtered.0.bed $SIZES ${OUTFILE%.bb}.all.bb 2>>$TMP/LOG.log
    echo >> $TMP/LOG.log

    (echo -n "Running bedToBigBed for filter A-2,reporting uniquelly mapped fragments after slippage filter: "; date) >> $TMP/LOG.log
    bedToBigBed $TMP/reads.filtered.2.bed $SIZES $OUTFILE 2>>$TMP/LOG.log

    if [ $UMI_flag = "TRUE" ]; then
    (echo -n "Running bedToBigBed for filter A-3, collapsed by UMI: "; date) >> $TMP/LOG.log
      bedToBigBed $TMP/reads.filtered.3.bed $SIZES ${OUTFILE%.bb}.UMI.bb 2>>$TMP/LOG.log
    #  rm $TMP/reads.filtered.3.bed
    fi


    #rm $TMP/reads.filtered.{0,2}.bed

else
    (echo -n "Running bedToBigBed: "; date) >> $TMP/LOG.log
    bedToBigBed $TMP/reads.filtered.bed $SIZES $OUTFILE 2>>$TMP/LOG.log
    #rm $TMP/reads.filtered.bed
fi

(echo -n "End: "
 date ) >> $TMP/LOG.log

if [ $KEEPLOG = "1" ]; then
    mv $TMP/LOG.log ${OUTFILE%.bb}.log
fi

if [ $KEEP = "1" ]; then
  mv $TMP/* ${OUTFILE%.bb}/.
fi

end=`date +%s`
runtime=$((end-start))
echo "done in $runtime"


################################################################################
# Exit
################################################################################

rm -rf ${TMP}

exit 0

