from Bio import SeqIO

import argparse
import re

# Regex para poder checar se sequencia eh ou nao uma sequencia normalizada
# ^atg(.*)(tga|taa|tag)$

# Para checar se a regex bate
# https://stackoverflow.com/questions/6576962/python-regular-expressions-return-true-false

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', dest="file")

def check_codons(fastafile):
    rgx_seqs = re.compile("^atg(.*)(tga|taa|tag)$")
    fasta_read = SeqIO.parse(fastafile, "fasta")
    for record in fasta_read:
        if rgx_seqs.match(str(record.seq)):
            print(record.format("fasta"))

args = parser.parse_args()
check_codons(args.file)
