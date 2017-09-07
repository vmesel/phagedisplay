import argparse
import os

from fasta_normalizer import *

# parser = argparse.ArgumentParser(description='Fasta line length normalizer')
#
# parser.add_argument('-r','--run',
#                     dest='run',
#                     action="store_true"
# )
#
# parser.add_argument('-f','--file',
#                     dest='file',
# )
#
# if __name__ == "__main__":
#     args = parser.parse_args()
#     if args.run:
#         fasta_one_liner(args.file, "schisto1.fasta")
#         os.system("sed -e 's/\*$//g' schisto1.fasta")

def fasta2dict(fasta_file):
    fasta_list = [line.strip() for line in open(fasta_file) if line.strip() != ""]
    fasta_dict = {}
    fasta_seq_names = [fasta_list[i].replace(">","") for i in range(0, len(fasta_list), 2)]
    fasta_sequences = [fasta_list[i] for i in range(1, len(fasta_list), 2)]

    for k, v in zip(fasta_seq_names, fasta_sequences):
        fasta_dict[k] = v

    return fasta_dict


fastadict = fasta2dict("../schisto2.fasta")

for k, v in fastadict.items():
    a = {"+": v.count("+"), "#": v.count("#")}
    if sum(a.values()) >= 2:
        print(">{}\n{}\n".format(k,v))