import argparse
import os

from fasta_normalizer import *

parser = argparse.ArgumentParser(description='Fasta line length normalizer')

parser.add_argument('-r','--run',
                    dest='run',
                    action="store_true"
)

parser.add_argument('-f','--file',
                    dest='file',
)

parser.add_argument('-s','--size',
                    dest='size',
)

if __name__ == "__main__":
    args = parser.parse_args()
    if args.run:
        wd = os.getcwd() + "/"
        fasta_one_liner(wd + args.file, "schisto1.fasta")
        os.system("sed -e 's/\*$//g' {wd}schisto1.fasta > {wd}schisto2.fasta".format(wd=wd))
        os.system("sed -e 's/:pep$//g' {wd}schisto2.fasta > {wd}schisto3.fasta".format(wd=wd))
        os.system("pyfasta split {wd}schisto3.fasta -n 1 -k {k} -o 7".format(wd=wd, k=args.size))
        seq_name_norm("{wd}schisto3.split.{s}mer.7overlap.fasta".format(wd=wd, s=args.size), "{wd}schisto4.fasta".format(wd=wd))
        fasta_normalizer(file="{wd}schisto4.fasta".format(wd=wd), desired_len=args.size, overlap_len=7, output_file = "{wd}schisto5.fasta".format(wd=wd))
        fasta_cleaner(input_fasta="{wd}schisto5.fasta".format(wd=wd), output_fasta= "{wd}schisto6.fasta".format(wd=wd), wrong_bases="X")




