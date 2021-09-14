
### Load arguments

import sys, getopt

### other parameters
# number of dinucleotide shuffled sequences per sequence as background
dinuc_shuffle_n=100

def main(argv):
   model_ID = ''
   try:
      opts, args = getopt.getopt(argv,"hm:s:c:",["model=", "sequence_set=", "class_output="])
   except getopt.GetoptError:
      print('DeepSTARR_nucl_contr_scores.py -m <CNN model file> -s <fasta seq file> -c <dev/hk>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('DeepSTARR_nucl_contr_scores.py -m <CNN model file> -s <fasta seq file> -c <dev/hk>')
         sys.exit()
      elif opt in ("-m", "--model"):
         model_ID = arg
      elif opt in ("-s", "--seq"):
         sequence_set = arg
      elif opt in ("-c", "--class_output"):
         class_output = arg
   if model_ID=='': sys.exit("CNN model file not found")
   if sequence_set=='': sys.exit("fasta seq file not found")
   if class_output=='': sys.exit("enhancer class (dev/hk) not found")
   print('CNN model file is ', model_ID)
   print('Input FASTA file is ', sequence_set)
   print('Enhancer class is ', class_output)
   return model_ID, sequence_set, class_output

if __name__ == "__main__":
   model_ID, sequence_set, class_output = main(sys.argv[1:])


### Load libraries
import pandas as pd

import sys
sys.path.append('Neural_Network_DNA_Demo/')
from helper import IOHelper, SequenceHelper # from https://github.com/bernardo-de-almeida/Neural_Network_DNA_Demo.git


### load fasta sequences functions
def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return

def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1


def prepare_input(file_seq):
    # Convert sequences to one-hot encoding matrix
    input_fasta_data_A = IOHelper.get_fastas_from_file(file_seq, uppercase=True)

    # get length of first sequence
    sequence_length = len(input_fasta_data_A.sequence.iloc[0])

    # Convert sequence to one hot encoding matrix
    seq_matrix_A = SequenceHelper.do_one_hot_encoding(input_fasta_data_A.sequence, sequence_length,
      SequenceHelper.parse_alpha_to_seq)
    print(seq_matrix_A.shape)
    
    X = np.nan_to_num(seq_matrix_A) # Replace NaN with zero and infinity with large finite numbers
    X_reshaped = X.reshape((X.shape[0], X.shape[1], X.shape[2]))
    
    print(file_seq)

    return X_reshaped

### load model functions
def load_model(model):
    import deeplift
    from keras.models import model_from_json
    keras_model_weights = model + '.h5'
    keras_model_json = model + '.json'
    keras_model = model_from_json(open(keras_model_json).read())
    keras_model.load_weights(keras_model_weights)
    return keras_model, keras_model_weights, keras_model_json

from deeplift.dinuc_shuffle import dinuc_shuffle
import numpy as np

### deepExplainer functions
def dinuc_shuffle_several_times(list_containing_input_modes_for_an_example,
                                seed=1234):
  assert len(list_containing_input_modes_for_an_example)==1
  onehot_seq = list_containing_input_modes_for_an_example[0]
  rng = np.random.RandomState(seed)
  to_return = np.array([dinuc_shuffle(onehot_seq, rng=rng) for i in range(dinuc_shuffle_n)])
  return [to_return] #wrap in list for compatibility with multiple modes

# get hypothetical scores also
def combine_mult_and_diffref(mult, orig_inp, bg_data):
    assert len(orig_inp)==1
    projected_hypothetical_contribs = np.zeros_like(bg_data[0]).astype("float")
    assert len(orig_inp[0].shape)==2
    #At each position in the input sequence, we iterate over the one-hot encoding
    # possibilities (eg: for genomic sequence, this is ACGT i.e.
    # 1000, 0100, 0010 and 0001) and compute the hypothetical 
    # difference-from-reference in each case. We then multiply the hypothetical
    # differences-from-reference with the multipliers to get the hypothetical contributions.
    #For each of the one-hot encoding possibilities,
    # the hypothetical contributions are then summed across the ACGT axis to estimate
    # the total hypothetical contribution of each position. This per-position hypothetical
    # contribution is then assigned ("projected") onto whichever base was present in the
    # hypothetical sequence.
    #The reason this is a fast estimate of what the importance scores *would* look
    # like if different bases were present in the underlying sequence is that
    # the multipliers are computed once using the original sequence, and are not
    # computed again for each hypothetical sequence.
    for i in range(orig_inp[0].shape[-1]):
        hypothetical_input = np.zeros_like(orig_inp[0]).astype("float")
        hypothetical_input[:,i] = 1.0
        hypothetical_difference_from_reference = (hypothetical_input[None,:,:]-bg_data[0])
        hypothetical_contribs = hypothetical_difference_from_reference*mult[0]
        projected_hypothetical_contribs[:,:,i] = np.sum(hypothetical_contribs,axis=-1) 
    return [np.mean(projected_hypothetical_contribs,axis=0)]


def my_deepExplainer(model, one_hot, class_output):
    import shap # forked from https://github.com/AvantiShri/shap/blob/master/shap/explainers/deep/deep_tf.py
    import numpy as np
    
    # output layer
    if class_output=="dev":
        out_layer=-2
    if class_output=="hk":
        out_layer=-1

    explainer = shap.DeepExplainer((model.layers[0].input, model.layers[out_layer].output),
        data=dinuc_shuffle_several_times,
        combine_mult_and_diffref=combine_mult_and_diffref)
        
    # running on all sequences
    shap_values_hypothetical = explainer.shap_values(one_hot)
    
    # normalising contribution scores
    # sum the deeplift importance scores across the ACGT axis (different nucleotides at the same position)
    # and “project” that summed importance onto whichever base is actually present at that position
    shap_values_contribution=shap_values_hypothetical[0]*one_hot

    return shap_values_hypothetical[0], shap_values_contribution


### Loading sequences and calculating nucleotide contribution scores

print("\nLoading sequences and model ...\n")

X_all = prepare_input(sequence_set)
keras_model, keras_model_weights, keras_model_json = load_model(model_ID)

print("\nRunning DeepExplain ...\n")

scores=my_deepExplainer(keras_model, X_all, class_output=class_output)

print("\nSaving ...\n")

import h5py
import os
model_ID_out=os.path.basename(model_ID)

if (os.path.isfile(sequence_set+"_"+model_ID_out+"_"+class_output+"_contribution_scores.h5")):
    os.remove(str(sequence_set+"_"+model_ID_out+"_"+class_output+"_contribution_scores.h5"))
f = h5py.File(sequence_set+"_"+model_ID_out+"_"+class_output+"_contribution_scores.h5")

g = f.create_group("contrib_scores")
# save the actual contribution scores
g.create_dataset(class_output, data=scores[1])
print("Done for " + class_output)

g = f.create_group("hyp_contrib_scores")
# save the hypothetical contribution scores
g.create_dataset(class_output, data=scores[0])
print("Done for " + class_output)

f.close()