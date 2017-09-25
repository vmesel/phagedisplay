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
        os.system("echo 'SCHISTO1' && cat {wd}schisto1.fasta | grep '>' | wc -l".format(wd=wd))
        os.system("sed -e 's/\*$//g' {wd}schisto1.fasta > {wd}schisto2.fasta".format(wd=wd))
        os.system("echo 'SCHISTO2' && cat {wd}schisto2.fasta | grep '>' | wc -l".format(wd=wd))
        os.system("sed -e 's/:pep$//g' {wd}schisto2.fasta > {wd}schisto3.fasta".format(wd=wd))
        os.system("echo 'SCHISTO3' && cat {wd}schisto3.fasta | grep '>' | wc -l".format(wd=wd))
        fasta_cleaner(input_fasta="{wd}schisto3.fasta".format(wd=wd), output_fasta= "{wd}schisto4.fasta".format(wd=wd), wrong_bases="X")
        os.system("echo 'SCHISTO4' && cat {wd}schisto4.fasta | grep '>' | wc -l".format(wd=wd))
        os.system("pyfasta split {wd}schisto4.fasta -n 1 -k {k} -o 7".format(wd=wd, k=args.size))
        os.system("echo 'SCHISTO4 splitted' && cat {wd}schisto4.split.{s}mer.7overlap.fasta | grep '>' | wc -l".format(wd=wd, s=args.size))
        seq_name_norm("{wd}schisto4.split.{s}mer.7overlap.fasta".format(wd=wd, s=args.size), "{wd}schisto5.fasta".format(wd=wd))
        os.system("echo 'SCHISTO5' && cat {wd}schisto5.fasta | grep '>' | wc -l".format(wd=wd))
        fasta_normalizer(file="{wd}schisto5.fasta".format(wd=wd), desired_len=args.size, overlap_len=7, output_file = "{wd}schisto6.fasta".format(wd=wd))
        os.system("echo 'SCHISTO6' && cat {wd}schisto6.fasta | grep '>' | wc -l".format(wd=wd))
