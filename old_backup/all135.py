from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', dest="file")

def complete_seqs(fasta_file):
    seqs = {}
 
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in seqs:
            seqs[record.id].append(record.description + "$$$" + str(record.seq))
        else:
            seqs[record.id] = []
            seqs[record.id].append(record.description + "$$$" + str(record.seq))

    ids_n_trabalhados = []
    for id in seqs:
        if len(seqs[id]) >= 2:
            incomp = seqs[id][-1].split("$$$")
            comp = seqs[id][-2].split("$$$")
            l = 135 - len(incomp[1])
            seq_comp = comp[1][-l:] + incomp[1]
            seqs[id][-1] = incomp[0] + "$$$" + seq_comp

    for id in seqs:
        for sequence in seqs[id]:
            seq = sequence.split("$$$")
            print(">" + seq[0].replace(" ", "_") + "\n" + seq[1] + "\n")

args = parser.parse_args()
complete_seqs(args.file)