# Scripts for training and interpreting DeepSTARR

## Training DeepSTARR

<p float="left" style="margin-bottom:0;margin-top:0;">
    <img height="200" src="../img/DeepSTARR.png">
    <img height="200" src="../img/DeepSTARR_predictions.png">
</p>

Code to train DeepSTARR is in the notebook [DeepSTARR_training](DeepSTARR_training.ipynb).
Data used to train and evaluate the DeepSTARR model as well as the final trained model are available on zenodo at https://doi.org/10.5281/zenodo.5502060.

### Predict developmental and housekeeping enhancer activity of new DNA sequences
To predict the developmental and housekeeping enhancer activity in *Drosophila melanogaster* S2 cells for new DNA sequences, please download the trained DeepSTARR model from [zenodo](https://doi.org/10.5281/zenodo.5502060) and run:
```
# Clone this repository
git clone git@github.com:bernardo-de-almeida/DeepSTARR.git
cd DeepSTARR

# create 'DeepSTARR' conda environment
conda env create -f DeepSTARR_conda_env.yml

# Activate the conda environment
source activate DeepSTARR

# Run prediction script
python DeepSTARR_pred_new_sequence.py -s Sequences_example.fa -m DeepSTARR.model
```
Where:
* -s FASTA file with input DNA sequences

## Interpreting DeepSTARR with nucleotide contribution scores
To compute nucleotide contribution scores for new DNA sequences in respect to developmental or housekeeping enhancer activity, please download the trained DeepSTARR model from [zenodo](https://doi.org/10.5281/zenodo.5502060) and run:
```
# Create and activate the conda environment as above, and then run:
python DeepSTARR_nucl_contr_scores.py -m DeepSTARR.model -s Sequences_example.fa -c dev # Contribution scores for developmental enhancer activity
python DeepSTARR_nucl_contr_scores.py -m DeepSTARR.model -s Sequences_example.fa -c hk # Contribution scores for housekeeping enhancer activity
```
Where:
* -s FASTA file with input DNA sequences
* -c Enhancer type for which contribution scores should be derived

#### Note
Neural_Network_DNA_Demo forked from https://github.com/const-ae/Neural_Network_DNA_Demo

## Questions
If you have any questions/requests/comments please contact me at [bernardo.almeida94@gmail.com](mailto:bernardo.almeida94@gmail.com).
