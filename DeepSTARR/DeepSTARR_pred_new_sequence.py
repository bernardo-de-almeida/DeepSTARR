
### Load arguments

import sys, getopt

def main(argv):
   new_seq = ''
   model_ID = ''
   try:
      opts, args = getopt.getopt(argv,"hs:m:",["seq=","model="])
   except getopt.GetoptError:
      print('DeepSTARR_pred_new_sequence.py -s <fasta seq file> -m <CNN model file>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('DeepSTARR_pred_new_sequence.py -s <fasta seq file> -m <CNN model file>')
         sys.exit()
      elif opt in ("-s", "--seq"):
         new_seq = arg
      elif opt in ("-m", "--model"):
         model_ID = arg
   if new_seq=='': sys.exit("fasta seq file not found")
   if model_ID=='': sys.exit("CNN model file not found")
   print('Input FASTA file is ', new_seq)
   print('Model file is ', model_ID)
   return new_seq, model_ID

if __name__ == "__main__":
   new_seq, model_ID = main(sys.argv[1:])


### Load libraries

from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers import BatchNormalization, InputLayer, Input
from keras.models import Sequential
from keras.optimizers import Adam
from keras.callbacks import EarlyStopping, History

import pandas as pd
import numpy as np

import sys
sys.path.append('Neural_Network_DNA_Demo/')
from helper import IOHelper, SequenceHelper # from https://github.com/bernardo-de-almeida/Neural_Network_DNA_Demo.git

### Load sequences
print("\nLoading sequences ...\n")
input_fasta = IOHelper.get_fastas_from_file(new_seq, uppercase=True)
print(input_fasta.shape)

# length of first sequence
sequence_length = len(input_fasta.sequence.iloc[0])

# Convert sequence to one hot encoding matrix
seq_matrix = SequenceHelper.do_one_hot_encoding(input_fasta.sequence, sequence_length,
                                                SequenceHelper.parse_alpha_to_seq)

### load model
def load_model(model_path):
    import deeplift
    from keras.models import model_from_json
    keras_model_weights = model_path + '.h5'
    keras_model_json = model_path + '.json'
    keras_model = model_from_json(open(keras_model_json).read())
    keras_model.load_weights(keras_model_weights)
    #keras_model.summary()
    return keras_model, keras_model_weights, keras_model_json

keras_model, keras_model_weights, keras_model_json = load_model(model_ID)

### predict dev and hk activity
print("\nPredicting ...\n")
pred=keras_model.predict(seq_matrix)
out_prediction = input_fasta
out_prediction['Predictions_dev'] = pred[0]
out_prediction['Predictions_hk'] = pred[1]

### save file
print("\nSaving file ...\n")
import os.path
model_ID_out=os.path.basename(model_ID)
out_prediction.to_csv(new_seq + "_predictions_" + model_ID_out + ".txt", sep="\t", index=False)
