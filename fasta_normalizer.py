from Bio import SeqIO


def fasta_normalizer(file, desired_len):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    fasta_dict_corrected = {}
    for k, v in fasta_dict.items():
        seq_name = v.description.replace(" ", "__")
        fasta_dict_corrected[seq_name] = []
        intermediate = str(v.seq)
        while True:
            fasta_dict_corrected[seq_name].append(intermediate[:desired_len])
            intermediate = intermediate[desired_len:]
            if intermediate == "":
                break
    for k, v in fasta_dict_corrected.items():
        if len(fasta_dict_corrected[k]) >= 2:
            complete_seq = fasta_dict_corrected[k][-2]
            incomplete_seq = fasta_dict_corrected[k][-1]
            left_to_complete = desired_len - len(incomplete_seq)
            completed_seq = complete_seq[-desired_len:] + incomplete_seq
            fasta_dict_corrected[k][-1] = completed_seq

    # for k, v in fasta_dict_corrected.items():


    return fasta_dict_corrected

a = fasta_normalizer("/Users/vmesel/GitHub/phagedisplay/Schisto.mRNA.Smp_and_product_name.2017.08.18.20.15.59.fasta", 135)
print(a)
