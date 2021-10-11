import numpy as np


def parse_alpha_to_seq(sequence):
    output = np.arange(len(sequence))
    for i in range(0, len(sequence)):
        snippet = sequence[i]
        if snippet == 'A':
            output[i] = 0
        elif snippet == 'C':
            output[i] = 1
        elif snippet == 'G':
            output[i] = 2
        elif snippet == 'T':
            output[i] = 3
        elif snippet == 'N':
            output[i] = -1
        else:
            raise AssertionError("Cannot handle snippet: " + snippet)
    return output


def parse_binary_seq(sequence):
    output = np.arange(len(sequence) / 2)
    for i in range(0, len(sequence), 2):
        snippet = sequence[i] + sequence[i + 1]
        if snippet == '00':
            output[int(i / 2)] = 0
        elif snippet == '01':
            output[int(i / 2)] = 1
        elif snippet == '10':
            output[int(i / 2)] = 2
        elif snippet == '11':
            output[int(i / 2)] = 3
        else:
            raise AssertionError("Cannot handle snippet: " + snippet)
    return output


def parse_binary_seq_to_alpha(sequence):
    output = ""
    for i in range(0, len(sequence), 2):
        snippet = sequence[i] + sequence[i + 1]
        if snippet == '00':
            output += 'A'
        elif snippet == '01':
            output += 'C'
        elif snippet == '10':
            output += 'G'
        elif snippet == '11':
            output += 'T'
        else:
            raise AssertionError("Cannot handle snippet: " + snippet)
    return output


def to_categorical(y, nb_classes=None):
    '''Convert class vector (integers from 0 to nb_classes)
    to binary class matrix, for use with categorical_crossentropy
    '''
    y = np.asarray(y, dtype='int32')
    if not nb_classes:
        nb_classes = np.max(y) + 1
    Y = np.zeros((len(y), nb_classes))
    for i in range(len(y)):
        if y[i] != -1:
            Y[i, y[i]] = 1.
    return Y


def do_one_hot_encoding(sequence, seq_length, f=parse_alpha_to_seq):
    X = np.zeros((sequence.shape[0], seq_length, 4))
    for idx in range(0, len(sequence)):
        X[idx] = to_categorical(f(sequence[idx]), 4)
    return X


def do_dinucleotide_shuffling(X, size=1):
    x_shuffled = np.repeat(X, size, 0)

    for x in range(0, x_shuffled.shape[0]):
        random_index = np.arange(0, X.shape[1]/2)
        np.random.shuffle(random_index)
        for y in range(0, int(X.shape[1]/2)):
            x_shuffled[x,y*2, ] = X[x%X.shape[0],random_index[y]*2]
            x_shuffled[x,(y*2)+1, ] = X[x%X.shape[0],(random_index[y]*2)+1]

    return x_shuffled


def generate_complementary_sequence(sequence):
    comp_seq = []
    for b in sequence:
        if b == "A":
            comp_seq.append("T")
        elif b == "T":
            comp_seq.append("A")
        elif b == "C":
            comp_seq.append("G")
        elif b == "G":
            comp_seq.append("C")
        elif b == "N":
            comp_seq.append("N")
        else:
            raise ValueError("Cannot convert base {0} to complement base!".format(b))
    return ''.join(comp_seq)