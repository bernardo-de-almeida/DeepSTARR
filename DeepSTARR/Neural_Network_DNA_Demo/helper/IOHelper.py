import gzip
import math
import os.path
from subprocess import Popen, PIPE, STDOUT

import numpy as np
import pandas as pd


def get_fastas_from_file(fasta_path, as_dict=False,
                         uppercase=False, stop_at=None):
    fastas = []
    seq = None
    header = None
    for r in (gzip.open(fasta_path) if fasta_path.endswith(".gz") else open(fasta_path)):
        if type(r) is bytes:
            r = r.decode("utf-8")
        r = r.strip()
        if r.startswith(">"):
            if seq != None and header != None:
                fastas.append([header, seq])
                if stop_at != None and len(fastas) >= stop_at:
                    break
            seq = ""
            header = r[1:]
        else:
            if seq != None:
                seq += r.upper() if uppercase else r
            else:
                seq = r.upper() if uppercase else r
    # append last fasta read by method
    if stop_at != None and len(fastas) < stop_at:
        fastas.append([header, seq])
    elif stop_at == None:
        fastas.append([header, seq])
    if as_dict:
        return {h: s for h, s in fastas}

    return pd.DataFrame({'location': [e[0] for e in fastas], 'sequence': [e[1] for e in fastas]})


def get_shape_fastas_from_file(fasta_path, as_dict=False,
                         uppercase=False, stop_at=None):
    fastas = []
    seq = None
    header = None
    for r in (gzip.open(fasta_path) if fasta_path.endswith(".gz") else open(fasta_path)):
        if type(r) is bytes:
            r = r.decode("utf-8")
        r = r.strip()
        if r.startswith(">"):
            if seq != None and header != None:
                fastas.append([header, seq])
                if stop_at != None and len(fastas) >= stop_at:
                    break
            seq = None
            header = r[1:]
        else:
            if seq != None:
                seq += "," + (r.upper() if uppercase else r)
            else:
                seq = r.upper() if uppercase else r
    # append last fasta read by method
    if stop_at != None and len(fastas) < stop_at:
        fastas.append([header, seq])
    elif stop_at == None:
        fastas.append([header, seq])
    if as_dict:
        return {h: s for h, s in fastas}

    return pd.DataFrame({'location': [e[0] for e in fastas], 'sequence': [e[1] for e in fastas]})


def get_padded_sequences(fasta_file):
    fasta = get_fastas_from_file(fasta_file)
    max_length = max([len(x) for x in fasta.sequence])
    padded_sequences = []
    for seq in fasta.sequence:
        diff = max_length - len(seq)
        n_seq = (math.floor(diff/2) * 'N') + seq + (math.ceil(diff/2) * 'N')
        padded_sequences.append(n_seq)
    fasta.sequence = padded_sequences
    return fasta


def convert_bed_to_fasta_hg19(bed_path, fasta_path, reference_genome_path, use_peak_max=False,
                              bp_flanking=50):
    '''
    Copied from Ignacio: /g/scb/zaugg/rio/EclipseProjects/zaugglab/lib/FastaAnalyzer.py
    :param bed_path: The path to our BED file
    :param fasta_path: The output fasta that will be created
    :param use_peak_max: If True, we will extract w.r.t. to peak position
    (See https://www.biostars.org/p/102710/ for format description
    :param bp_flanking: If use_peak is True, then flanking regions will
    be calculated from this file
    :return:
    '''

    args = ["/g/software/bin/bedtools", "getfasta", "-fi", reference_genome_path,
            "-fo", fasta_path]

    # create a new coordinates file with flanking sequences
    if use_peak_max:
        df = pd.read_csv(bed_path, sep='\t', index_col=False,
                         names=['chrom', 'chromStart', 'chromEnd', 'name', 'score',
                                'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
        df['startFromPeak'] = df['thickStart'] - bp_flanking
        df['endFromPeak'] = df['thickStart'] + bp_flanking
        df = df[['chrom', 'startFromPeak', 'endFromPeak']]
        tsv_string = df.to_csv(header=False, sep='\t', index=False)
        args = args + ['-bed', 'stdin']

        p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        x = p.communicate(input=tsv_string.encode(encoding='UTF-8'))
        x = x[0].decode('UTF-8')
        if x != '':
            print("ERROR: " + x)
    else:
        os.system(" ".join(args + ['-bed', bed_path]))


def write_fasta_file(file, sequences, descr=None):
    """
    Sequences has to be a list of strings. descr can be None than a dummy line is inserted or a list of the
    same length as sequences.
    """
    with open(file, "w") as out:
        for idx, seq in enumerate(sequences):
            if descr is None:
                out.write(">Dummy_Line\n")
            else:
                out.write(">" + str(descr[idx]) + "\n")
            out.write("".join(seq) + "\n")


def save_keras_model(model, model_path, overwrite=False):
    json_string = model.to_json()
    with open(model_path + '.json', 'w+') as f:
        f.write(json_string)
    model.save_weights(model_path + '.h5', overwrite=overwrite)


def load_keras_model(path):
    from keras.models import model_from_json
    model = model_from_json(open(path + '.json').read())
    model.load_weights(path + '.h5')
    return model


def save_scoring_file(header, values, scores, labels, file):
    if len(scores) != len(labels):
        raise ValueError("The score and label length must match!")
    if len(header) != scores.shape[3] + values.shape[2]:
        raise ValueError("The value + score width and header length must match!")

    with open(file, 'w') as output:
        output.write("\t".join(["Index", "Label"] + header) + "\n")
        for line_idx in range(0, len(scores)):
            output.write("\t".join([str(line_idx), labels[line_idx]]
                                   + ["["+",".join(map(str, values[line_idx, :, c]))+"]" for c in range(0, values.shape[2])]
                                   + ["["+",".join(map(str, scores[line_idx, 0, :, c]))+"]" for c in range(0, scores.shape[3])]))
            output.write("\n")



def read_importance_file(location):
    return pd.read_csv(location, sep="\t")


def parse_importance_df(df, col_names):
    # Iterate over every entry
    parsed_cols = []
    for name in col_names:
        col = df[name].as_matrix()
        parsed_col = np.apply_along_axis(lambda e: np.array([float(x) for x in e[0][1:-1].split(",")]), 1, col.reshape(len(col),1))
        parsed_cols.append(parsed_col)
    return np.stack(parsed_cols, 2)


def write_output_file(output_file, name, PositiveData, NegativeData, Training_Script, aucs, auprcs, importance_scores):
    with open(output_file, "w") as out:
        out.write("Name:" + str(name) + "\n")
        out.write("PositiveData:" + str(PositiveData) + "\n")
        out.write("NegativeData:" + str(NegativeData) + "\n")
        out.write("Training_Script:" + str(Training_Script) + "\n")
        out.write("AUCs:" + ",".join(map(str, aucs)) + "\n")
        out.write("AUPRCs:" + ",".join(map(str, auprcs)) + "\n")
        out.write("Importance_Scores:" + str(importance_scores) + "\n")