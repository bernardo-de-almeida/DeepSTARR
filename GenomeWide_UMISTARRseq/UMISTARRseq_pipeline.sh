#######################
## Process genome-wide STARR-seq sequencing data
#######################

folder=GenomeWide_UMISTARRseq
cd $folder

# folder to write results
dataFolder=$folder/data
mkdir -p $dataFolder

# wrapper to submit jobs to cluster
bsub=bsub_gridengine

## NOTE ##
# The original sequencing data at our institute is provided in a BAM file containing the reads from a whole lane together with their the i5 and i7 indexes.
# The script below is based on this BAM file as input and requires barcodes to demultiplex the reads and map them with bowtie.
# However we can only provide fastq files of the demultiplexed reads for each experiment.
# For samples with unique molecular identifiers (UMIs) at the i7 index, the UMI information is included in the read name and can be used to collapse reads with identical UMIs.

#######################
## paired-end mapping with bowtie
#######################

## make an exeriment file with all sample information and respective barcodes for demultiplexing from main sequencing BAM file --> $dataFolder/experiment.txt

head $dataFolder/experiment.txt
experimentFile=$dataFolder/experiment.txt

GENOME="dm3" #download genome fasta file (no_chrU_chrUextra_chrM)
mkdir -p log_mapping

ALIGNER=bowtie_pe.sh

grep -v -E "^#" $experimentFile | \
    while read line; do
    INFILE=$( echo $line | awk '{print $2}' )
    outdir=$dataFolder
    BARCODES14=$( echo $line | awk '{print $9}' )
    BARCODES12=$( echo $line | awk '{print $11}' )
    OUTFILE=$( echo $line | awk '{print $10}')
    BARCODE14_LEN=$( echo $BARCODES14 | awk '{print length($1)}' )
    BARCODE12_LEN=$( echo $BARCODES12 | awk '{print length($1)}' )
    
    # with UMIs
    if [ "$BARCODES14" == "UMI" ]; then
      BARCODE14_LEN=10 # length of UMI
      $bsub -o log_mapping -C 10 -T '2-00:00:00' -n "${OUTFILE}_mapping" "${ALIGNER} -i $INFILE -o ${outdir}/${OUTFILE}.bb -B $BARCODES12 -L $BARCODE12_LEN -l $BARCODE14_LEN --umi -f A -g $GENOME " > log_mapping/msg.$OUTFILE.tmp
    else
      # no UMIs
      $bsub -o log_mapping -C 10 -n "${OUTFILE}_mapping" "${ALIGNER} -i $INFILE -o ${outdir}/${OUTFILE}.bb -B $BARCODES12 -L $BARCODE12_LEN -b $BARCODES14 -l $BARCODE14_LEN -f A -g $GENOME " > log_mapping/msg.$OUTFILE.tmp
    fi

    # get job IDs to wait for mapping to finish
    ID=$(paste log_mapping/msg.$OUTFILE.tmp|grep Submitted|awk '{print $4}')
    if [ "$OS_id" == "Debian" ]; then
      ID="${OUTFILE}_mapping"
    fi

    # bigwig files
    $bsub -o log_bigwig -C 5 -d "$ID" "bigBedToBigWig.sh -i ${outdir}/${OUTFILE}.all.bb -t libSize -g $GENOME" #all mapped reads
    $bsub -o log_bigwig -C 5 -d "$ID" "bigBedToBigWig.sh -i ${outdir}/${OUTFILE}.UMI.bb -t libSize -g $GENOME" #UMI-collapsed reads

done

#######################
## select reads with specific lengths - 150-250bp
#######################

# bigBedToBed from kent_tools

SIZES=dm3.chrom.sizes

for f in `ls ${dataFolder}/*200*bb`
  do echo $f
  bigBedToBed $f stdout | awk '{if($3-$2 >= 150 && $3-$2 <= 250)print $0}' > test.bed
  bedToBigBed test.bed $SIZES ${f%.bb}_cut.bb
  ${prog_dir}/bigBedToBigWig.sh -i ${f%.bb}_cut.bb -t libSize -g $GENOME
  rm test.bed
done

#######################
## call STARR-seq peaks
#######################

### join replicates to call peaks
for f in DSCP_200bp_gw RpS12_200bp_gw
  do echo $f
  bigBedToBed ${dataFolder}/${f}_Rep1.UMI_cut.bb rep1
  bigBedToBed ${dataFolder}/${f}_Rep2.UMI_cut.bb rep2
  cat rep1 rep2 | sort -k1,1 -k2,2n > sorted.bed
  bedToBigBed sorted.bed $SIZES ${dataFolder}/${f}.UMI_cut_merged.bb
  ${prog_dir}/bigBedToBigWig.sh -i ${dataFolder}/${f}.UMI_cut_merged.bb -t libSize -g $GENOME
  rm rep1
  rm rep2
  rm sorted.bed
done

### Get median fragment size
rm ${dataFolder}/fragment_size_median.txt
for f in `ls ${dataFolder}/*_cut*.bb`
  do echo $f
  size=$(bigBedToBed $f stdout | awk '{ print $3-$2+1 }' | median.R -i - )
  echo -e "$(basename $f)\t$size" >> ${dataFolder}/fragment_size_median.txt
done

### call peaks
for f in DSCP_200bp_gw.UMI RpS12_200bp_gw.UMI
  do echo $f
  if [ "$f" == "DSCP_200bp_gw.UMI" ]; then Input=input_DSCP_200bp.all_cut.bb; fi
  if [ "$f" == "RpS12_200bp_gw.UMI" ]; then Input=input_RPS12_200bp.all_cut.bb; fi
  echo $Input
  bsub -o log_peaks -C 5 "module load python/2.7.13-foss-2017a; module load rpy2/2.8.6-foss-2017a-python-2.7.13; \
  call_peaks.sh -e ${dataFolder}/${f}_cut_merged.bb -b ${dataFolder}/$Input -g ${GENOME} -M 3 -W 1000 -Z 1.67 -P 0.001"
done
