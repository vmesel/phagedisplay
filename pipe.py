import argparse
import os

from fasta_normalizer import *
from codon_translation import *
from motif_remover import *

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

parser.add_argument('-o','--overlap',
                    dest='overlap',
)

if __name__ == "__main__":
    args = parser.parse_args()
    if args.run:
        wd = os.getcwd() + "/"
        fasta_one_liner(wd + args.file, "schisto1.fasta")
        os.system("echo 'SCHISTO1' && cat {wd}schisto1.fasta | grep '>' | wc -l".format(wd=wd))
        os.system("sed -e 's/\*$//g' {wd}schisto1.fasta > {wd}schisto2.fasta".format(wd=wd)) # REMOVE STOP CODONS
        os.system("echo 'SCHISTO2' && cat {wd}schisto2.fasta | grep '>' | wc -l".format(wd=wd))
        os.system("sed -e 's/:pep$//g' {wd}schisto2.fasta > {wd}schisto3.fasta".format(wd=wd)) # LIMPA NOME DA SEQUENCIA
        os.system("echo 'SCHISTO3' && cat {wd}schisto3.fasta | grep '>' | wc -l".format(wd=wd))
        fasta_cleaner(input_fasta="{wd}schisto3.fasta".format(wd=wd), output_fasta= "{wd}schisto4.fasta".format(wd=wd), wrong_bases="X") # LIMPA A SEQUENCIA
        os.system("echo 'SCHISTO4' && cat {wd}schisto4.fasta | grep '>' | wc -l".format(wd=wd))
        # AQUI DEVE IR O BLAST
        # AQUI DEVE VIR A REMOCAO DE DUPLICATAS DO BLAST
        # AQUI DEVE VIR A REMOCAO DE 906 PROBES
        fasta_translation_ecoli(fasta_file="{}schisto4.fasta".format(wd), output_file="{}schisto5.fasta".format(wd)) # TRADUZ PARA E. COLI
        os.system("pyfasta split {wd}schisto5.fasta -n 1 -k {k} -o {o}".format(wd=wd, k=args.size, o=args.overlap))
        os.system("echo 'SCHISTO5 splitted' && cat {wd}schisto5.split.{s}mer.{o}overlap.fasta | grep '>' | wc -l".format(wd=wd, s=args.size, o=args.overlap))
        seq_name_norm("{wd}schisto5.split.{s}mer.{o}overlap.fasta".format(wd=wd, s=args.size, o=args.overlap), "{wd}schisto6.fasta".format(wd=wd))
        os.system("echo 'SCHISTO6' && cat {wd}schisto6.fasta | grep '>' | wc -l".format(wd=wd))
        fasta_normalizer(file="{wd}schisto6.fasta".format(wd=wd), desired_len=args.size, overlap_len=args.overlap, output_file = "{wd}schisto7.fasta".format(wd=wd))
        os.system("echo 'SCHISTO7' && cat {wd}schisto7.fasta | grep '>' | wc -l".format(wd=wd))
        call_cross_match("{wd}schisto7.fasta".format(wd=wd), "nco1", "{wd}schisto8_nco1.fasta".format(wd=wd))
        print("SCHISTO8 - NCO1")
        call_cross_match("{wd}schisto8_nco1.fasta".format(wd=wd), "pst1", "{wd}schisto8_nco1_pst1.fasta".format(wd=wd))
        print("SCHISTO8 - NCO1 COM PST1")
        call_cross_match("{wd}schisto8_nco1_pst1.fasta".format(wd=wd), "ecor5", "{wd}schisto8_nco1_pst1_ecor5.fasta".format(wd=wd))
        print("SCHISTO8 - NCO1, PST1 COM ECOR5")

        # ATE AQUI EH COM O CROSS MATCH
        call_cross_match("{wd}schisto8_nco1_pst1_ecor5.fasta".format(wd=wd), "5a", "{wd}schisto8_nco1_pst1_ecor5_5a.fasta".format(wd=wd))
        print("SCHISTO8 - NCO1, PST1, ECOR5 COM 5A")
        call_cross_match("{wd}schisto8_nco1_pst1_ecor5_5a.fasta".format(wd=wd), "5t", "{wd}schisto8_wout_motifs.fasta".format(wd=wd))
        print("SCHISTO8 - NCO1, PST1, ECOR5, 5A E 5T")
        call_cross_match("{wd}schisto8_wout_motifs.fasta".format(wd=wd), "nco1", "{wd}schisto9_nco1.fasta".format(wd=wd))
        print("SCHISTO9 - 5A")
        call_cross_match("{wd}schisto9_nco1.fasta".format(wd=wd), "pst1", "{wd}schisto9_nco1_pst1.fasta".format(wd=wd))
        print("SCHISTO9 - 5A")
        call_cross_match("{wd}schisto9_nco1_pst1.fasta".format(wd=wd), "ecor5", "{wd}schisto9_nco1_pst1_ecor5.fasta".format(wd=wd))
        print("SCHISTO9 - 5A")
        call_cross_match("{wd}schisto9_nco1_pst1_ecor5.fasta".format(wd=wd), "5a", "{wd}schisto9_5a.fasta".format(wd=wd))
        print("SCHISTO9 - 5A")
        call_cross_match("{wd}schisto9_5a.fasta".format(wd=wd), "5t", "{wd}schisto9_wout_motifs.fasta".format(wd=wd))
        print("SCHISTO9 - 5A E 5T")


        # AQUI EH A PARTE COM O REPLACE

        # FAZER PARTE DOS ADAPTADORES
