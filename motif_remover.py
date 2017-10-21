import itertools
from fasta_normalizer import fasta2dict, dict2fasta, fasta_one_liner
import os

def frame_calc(sequence):
    modulo = len(sequence) % 3
    if modulo == 0:
        return 0
    elif modulo == 1:
        return 1
    elif modulo == 2:
        return 2


def motif_frame_counter(sequence):
    intermediario = []
    output = []
    for n, i in enumerate(list(sequence)):
        if i == "x" and n != 1:
            intermediario.append(n)
        elif i == "x" and n == 1:
            intermediario.append(0)

        if len(intermediario) == 6:
            output.append(frame_calc(sequence[:intermediario[0]]))
            intermediario = []

    return output


def replacer(x, rep0, rep1, rep2):
    if x == 0:
        return rep0
    elif x == 1:
        return rep1
    elif x == 2:
        return rep2


def motif_replacer(seq, rep0, rep1, rep2):
    sequence_wout_RS = seq.split("xxxxxx")
    sequence_to_replace = [
            replacer(x, rep0, rep1, rep2) for x in motif_frame_counter(seq)
    ]
    zip_sequence_out = list(itertools.zip_longest(sequence_wout_RS, sequence_to_replace, fillvalue=""))
    sequence_out = "".join([x[0] + x[1] for x in zip_sequence_out])
    return sequence_out


def fasta_motif_replacer(fasta, restriction_site, output_file):
    fasta_to_dict = fasta2dict(fasta, pyfastalib=False)
    rs_replacing = {
            "nco1":("ccttgg", "ctatgg", "ccacgg"),
            "pst1":("ctgcaa", "ctgcgg", "ctgtag"),
            "ecor5":("gatatt", "gatacc", "##ERROR1##"),
            "5a": ("aagaa", "aaaag", "aaaga"),
            "5t": ("ttctt", "ttttc", "tttct"),
    }
    rep0, rep1, rep2 = rs_replacing[restriction_site]
    dict_to_fasta = {k: motif_replacer(v.lower(), rep0, rep1, rep2).upper() for k, v in fasta_to_dict.items()}
    dict2fasta(dict_to_fasta, output_fasta=output_file)

def call_cross_match(seq_file, restriction_site, output_fasta):
    rs_replacing = {
        "nco1": "/work/users/vinicius/def_pipe/phagedisplay/motifs/nco1.vec",
        "pst1": "/work/users/vinicius/def_pipe/phagedisplay/motifs/pst1.vec",
        "ecor5": "/work/users/vinicius/def_pipe/phagedisplay/motifs/ecor5.vec",
        "5a": "/work/users/vinicius/def_pipe/phagedisplay/motifs/5a.vec",
        "5t": "/work/users/vinicius/def_pipe/phagedisplay/motifs/5t.vec",
    }
    motif_file = rs_replacing[restriction_site]
    if restriction_site not in ["5a","5t"]:
        motif_runner = "/home/elton/bioinformatics-tools/phredPhrapConsed/phrap/cross_match.manyreads {sf} {mf}  -minmatch 6 -minscore 6 -screen".format(sf=seq_file, mf=motif_file)
        os.system(motif_runner)
        to_one_liner = (seq_file + ".screen", seq_file + ".one_liner")
        fasta_one_liner(to_one_liner[0], to_one_liner[1])
        fasta_motif_replacer(to_one_liner[1], restriction_site, output_fasta)
    else:
        of = seq_file + ".pre_screen"
        # if restriction_site == "5a":
        #     os.system("sed 's/AAAAA/XXXXX/g' {inf} > {of}".format(inf=seq_file, of=of))
        # elif restriction_site == "5t":
        #     os.system("sed 's/TTTTT/XXXXX/g' {inf} > {of}".format(inf=seq_file, of=of))
        motif_runner = "/home/elton/bioinformatics-tools/phredPhrapConsed/phrap/cross_match.manyreads {sf} {mf}  -minmatch 5 -screen".format(sf=of, mf=motif_file)
        to_one_liner = (of + ".screen", of + ".one_liner")
        os.system(motif_runner)
        fasta_one_liner(to_one_liner[0], to_one_liner[1])
        fasta_motif_replacer(to_one_liner[1], restriction_site, output_fasta)
