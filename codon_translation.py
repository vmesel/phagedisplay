from Bio import SeqIO
from random import choice

def reverse_translate(sequence):
    translation_table = {
        "F":["ttt","ttc",],
        "L":["tta","ttg","ctt","ctc","ctg",],
        "I":["att","atc",],
        "M":["atg"],
        "V":["gtt","gtc","gta","gtg",],
        "S":["tct","tca","tcg","agt","agc",],
        "P":["cct","cca","ccg",],
        "T":["act","acc","acg",],
        "A":["gct","gcc","gca","gcg",],
        "Y":["tat","tac"],
        "H":["cat","cac",],
        "Q":["caa","cag",],
        "N":["aat","aac",],
        "K":["aag","aaa",],
        "D":["gat","gac",],
        "E":["gaa","gag",],
        "C":["tgt","tgc",],
        "W":["tgg"],
        "R":["cgt","cgc",],
        "G":["ggt","ggc","gga",],
	"U":[""],
    }
 
    final_sequence = []
    for nt in list(sequence):
        final_sequence.append(choice(translation_table[nt]))

    return "".join(final_sequence)



def fasta_translation_ecoli(fasta_file, output_file):
    w = open(output_file, "w+")
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    converted_fasta = {}
    for record, v in record_dict.items():
        converted_fasta[record] = reverse_translate(str(v.seq))

    for k, v in converted_fasta.items():
        w.write(">{k}\n{v}\n".format(k=k,v=v))

    w.close()
