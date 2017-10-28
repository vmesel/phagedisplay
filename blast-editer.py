import argparse
import os

import pandas as pd
from pprint import pprint

parser = argparse.ArgumentParser(description="BLAST File extract unique sequences")

parser.add_argument('-f', '--file', dest="file")

args = parser.parse_args()

wd = os.getcwd() + "/"

df = pd.read_csv(wd + args.file, sep="\t", header=None)
# print(df)
dict_smps = {}

for row in df.itertuples():
    if row[1] in dict_smps:
        if row[2] not in dict_smps[row[1]]:
            dict_smps[row[1]].append(row[2])
    else:
        dict_smps[row[1]] = []
        if row[2] not in dict_smps[row[1]]:
            dict_smps[row[1]].append(row[2])

valores_dict = []
list_to_unfold = []

for k, v in dict_smps.items():
    if v not in valores_dict:
        valores_dict.append(v)
        list_to_unfold.append((k,v))

for k,v in list_to_unfold:
    print("{k}\t{v}".format(k=k,v=v))