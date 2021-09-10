
#######################
## Process oligo genome-wide STARR-seq sequencing data
#######################

folder=Oligo_UMISTARRseq
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

# to submit to the cluster
mkdir -p log_mapping

# GENOME is a reference index containing 249 bp long sequences included in the library
ALIGNER=bowtie_pe_NoRevCompMapping.sh # no mapping against the reverse-complement reference strand

# map with no mismatches
BOWTIE_MM=0

grep -v -E "^#" $experimentFile | \
    while read line; do
    INFILE=$( echo $line | awk '{print $2}' )
    outdir=$dataFolder
    BARCODES14=$( echo $line | awk '{print $9}' )
    BARCODES12=$( echo $line | awk '{print $11}' )
    OUTFILE=$( echo $line | awk '{print $10}')
    BARCODE14_LEN=$( echo $BARCODES14 | awk '{print length($1)}' )
    BARCODE12_LEN=$( echo $BARCODES12 | awk '{print length($1)}' )
    GENOME=$( echo $line | awk '{print $5}' ) # using my own index

    # with UMIs
    if [ "$BARCODES14" == "UMI" ]; then
      BARCODE14_LEN=10 # length of UMI
      $bsub -o log_mapping -C 10 -T '5:00:00' -n "${OUTFILE}_mapping" "$ALIGNER -i $INFILE -o ${outdir}/${OUTFILE}.bb -B $BARCODES12 -L $BARCODE12_LEN -l $BARCODE14_LEN --umi -f A -g $GENOME -m $BOWTIE_MM" > log_mapping/msg.$OUTFILE.tmp
    else
      # no UMIs
      $bsub -o log_mapping -C 10 -n "${OUTFILE}_mapping" "$ALIGNER -i $INFILE -o ${outdir}/${OUTFILE}.bb -B $BARCODES12 -L $BARCODE12_LEN -b $BARCODES14 -l $BARCODE14_LEN -f A -g $GENOME -m $BOWTIE_MM" > log_mapping/msg.$OUTFILE.tmp
    fi
done
