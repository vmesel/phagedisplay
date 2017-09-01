import argparse

from Bio import SeqIO, Seq
from Bio.Alphabet import generic_dna
from Bio.Data.CodonTable import TranslationError

parser = argparse.ArgumentParser(description='Fasta line length normalizer')

parser.add_argument('-n','--fasta_normalizer',
                    dest='fasta_normalizer',
                    action="store_true"
)
parser.add_argument('-p','--peptide-fasta',
                    dest='peptide_fasta',
                    action="store_true"
)
parser.add_argument('-c','--fasta-cleaner',
                    dest='fasta_cleaner',
                    action="store_true"
)
parser.add_argument('-f','--file',
                    dest='fasta_file'
)
parser.add_argument('-l','--length',
                    dest='fasta_line_length'
)
parser.add_argument('-o','--overlap',
                    dest='fasta_sequences_overlap'
)
parser.add_argument('-of','--output_file',
                    dest='output_file'
)

def fasta2dict(fasta_file):
    fasta_list = [line.strip() for line in open(fasta_file) if line.strip() != ""]
    fasta_dict = {}
    fasta_seq_names = [fasta_list[i].replace(">","") for i in range(0, len(fasta_list), 2)]
    fasta_sequences = [fasta_list[i] for i in range(1, len(fasta_list), 2)]

    for k, v in zip(fasta_seq_names, fasta_sequences):
        splited_seq_name = k.split("|")
        seq_name = splited_seq_name[0]
        seq_pos = splited_seq_name[1]
        if seq_name in fasta_dict:
            fasta_dict[seq_name].append((seq_pos, v))
        else:
            fasta_dict[seq_name] = []
            fasta_dict[seq_name].append((seq_pos, v))

    return fasta_dict

def dict2fasta(fasta_cleaned, output_fasta=None):
    fasta_list = []
    for k, v in fasta_cleaned.items():
        for ind, seq in v:
            fasta_list.append(">{}|{}\n{}\n".format(k, ind, seq))

    if output_fasta != None:
        with open(output_fasta, "w+") as f:
            f.write("\n".join(fasta_list))
    else:
        print("\n".join(fasta_list))

def fasta_normalizer(file, desired_len, overlap_len, output_file = None):
    fasta_dict = fasta2dict(file)

    for k, v in fasta_dict.items():
        if len(v) > 1:
            seq_complete = v[-2][1]
            seq_incomplete = v[-1][1]

            if len(seq_incomplete) >= int(overlap_len):
                if len(seq_incomplete) < int(desired_len):
                    missing_len = int(desired_len) - len(seq_incomplete)
                    v[-1] = (v[-1][0], seq_complete[missing_len:] + seq_complete)

    fasta_file_output = []
    for k, v in fasta_dict.items():
        for ind, seq in v:
            if len(seq) == int(desired_len):
                fasta_file_output.append(">{}|{}\n{}\n".format(k, ind, seq))

    if output_file != None:
        with open(output_file, "w+") as f:
            f.write("\n".join(fasta_file_output))
    else:
        print("\n".join(fasta_file_output))


def fasta_cleaner(input_fasta, output_fasta=None):
    fasta_dict = fasta2dict(input_fasta)
    fasta_cleaned = {}

    for k, v in fasta_dict.items():
        for n, item in enumerate(v):
            if "n" not in item[-1]:
                if k in fasta_cleaned:
                    fasta_cleaned[k].append(item)
                else:
                    fasta_cleaned[k] = []
                    fasta_cleaned[k].append(item)

    dict2fasta(fasta_cleaned)

def pep_fasta(input_fasta, output_fasta=None):
    fasta_dict = fasta2dict(input_fasta)
    fasta_dict_bkp = {}
    for k, v in fasta_dict.items():
        # import ipdb; ipdb.set_trace()
        try:
            fasta_dict_bkp[k] = [(ind, Seq.Seq(x, generic_dna).translate()) for ind, x in v]
        except:
            pass

    dict2fasta(fasta_dict_bkp, output_fasta)

if __name__ == '__main__':
    args = parser.parse_args()
    if args.fasta_normalizer:
        fasta_normalizer(args.fasta_file,
        args.fasta_line_length,
        args.fasta_sequences_overlap,
        args.output_file
)

    if args.fasta_cleaner:
        fasta_cleaner(args.fasta_file,
            args.output_file
)

    if args.peptide_fasta:
        pep_fasta(args.fasta_file,
            args.output_file
)
