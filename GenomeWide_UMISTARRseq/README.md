# Scripts for processing genome-wide UMI-STARR-seq data

Main script: UMISTARRseq_pipeline.sh

Steps:
- Map reads with bowtie and collapse by UMIs
- Select reads with specific lengths - 150-250bp
- Create coverage BigWig files
- Call STARR-seq peaks

The raw sequencing and processed data are available from GEO under accession number GSE183939.
Genome browser tracks are available at https://genome.ucsc.edu/s/bernardo.almeida/DeepSTARR_manuscript.
