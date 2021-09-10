# DeepSTARR
DeepSTARR is a deep learning model built to quantitatively predict the activities of developmental and housekeeping enhancers from DNA sequence in *Drosophila melanogaster* S2 cells.

For more information, see the DeepSTARR manuscript:  
*DeepSTARR predicts enhancer activity from DNA sequence and enables the de novo design of enhancers* bioRxiv (2021)

This repository contains the code used to process genome-wide and oligo UMI-STARR-seq and train DeepSTARR.
The raw sequencing data are available from GEO under accession number GSE183939.

## Genome-wide enhancer activity maps of developmental and housekeeping enhancers
We used UMI-STARR-seq ([Arnold et al., 2013](http://www.sciencemag.org/lookup/doi/10.1126/science.1232542); [Neumayr et al., 2019](https://doi.org/10.1002/cpmb.105)) to generate genome-wide high resolution, quantitative activity maps of developmental and housekeeping enhancers, representing the two main transcriptional programs in *Drosophila* S2 cells ([Arnold et al., 2017](http://dx.doi.org/doi:10.1038/nbt.3739); [Haberle et al., 2019](https://doi.org/10.1038/s41586-019-1210-7); [Zabidi et al., 2015](http://dx.doi.org/10.1038/nature13994)).

<img src="img/gw_UMISTARRseq.png" width="700" style="margin-bottom:0;margin-top:0;"/>

You can find the code to process the data [here](GenomeWide_UMISTARRseq)

## DeepSTARR model


<p float="left" style="margin-bottom:0;margin-top:0;">
    <img height="200" src="img/DeepSTARR.png">
    <img height="200" src="img/DeepSTARR_predictions.png">
</p>

You can find the code [here](DeepSTARR)

### Predict developmental and housekeeping activity of new DNA sequences



## UMI-STARR-seq with designed oligo libraries to test more than 40,000 wildtype and mutant Drosophila and human enhancers

You can find the code [here](Oligo_UMISTARRseq)




## UCSC Genome Browser tracks
Genome browser tracks showing genome-wide UMI-STARR-seq and DeepSTARR predictions in Drosophila, together with the enhancers used for mutagenesis, mutated motif instances and respective log2 fold-changes, are available at https://genome.ucsc.edu/s/bernardo.almeida/DeepSTARR_manuscript.

## Questions
If you have any questions/requests/comments please contact me at [bernardo.almeida94@gmail.com](mailto:bernardo.almeida94@gmail.com)
