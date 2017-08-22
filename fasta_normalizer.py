from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Fasta line length normalizer')

parser.add_argument('-f','--file', dest='fasta_file')
parser.add_argument('-l','--length', dest='fasta_line_length')

def fasta_normalizer(file, desired_len):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    fasta_dict_corrected = {}
    for k, v in fasta_dict.items():
        seq_name = v.description.replace(" ", "__")
        fasta_dict_corrected[seq_name] = []
        intermediate = str(v.seq)
        while True:
            fasta_dict_corrected[seq_name].append(intermediate[:int(desired_len)])
            intermediate = intermediate[int(desired_len):]
            if intermediate == "":
                break

    for k, v in fasta_dict_corrected.items():
        print ">" + k
        print "\n".join(v) + "\n"


args = parser.parse_args()
fasta_normalizer(args.fasta_file, args.fasta_line_length)
