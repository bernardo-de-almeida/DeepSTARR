# Scripts for processing oligo UMI-STARR-seq data

Pipeline for mapping reads with bowtie (no mismtaches) to reference index (oligos included in the library) and UMI collapsing: [oligo_UMISTARRseq_pipeline.sh](oligo_UMISTARRseq_pipeline.sh)

R markdowns for processing mapped reads, check quality of screens, replicates, and calculate activity of each oligo with DESeq2:
- Drosophila library: [Drosophila_oligo_library_processing.Rmd](Drosophila_oligo_library_processing.Rmd)
- Human library: [Human_oligo_library_processing.Rmd](Human_oligo_library_processing.Rmd)

The raw sequencing and processed data are available from GEO under accession number [GSE183939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183939).
