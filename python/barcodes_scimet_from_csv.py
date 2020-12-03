## this script takes the csv file containing the barcode sequences and make a .txt file compatible for use with sciMET demultiplex script (perl)
## arguments needed: csv file containing the index name and sequence; output location
## the sequences should be the reverse complement of the sequences in the csv file

import sys
import pandas as pd

## import csv file
csv_file = sys.argv[1]
output_location = sys.argv[2]

pd.set_option("display.max_rows", 400) # this sets up the number of rows displayed
pd.set_option("display.float_format", lambda x: '%.3f' % x)
df = pd.read_csv(csv_file)
df = df[["Name", "Index"]]
df = df[df['Index'].notna()]

def rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

## create name:sequence dictionary
tn5_dict = {}
i5_dict = {}
i7_dict = {}

for row, index in df.iterrows():
    if index[0].startswith("sciMET_T"):
        tn5_dict[index[0]] = rev_complement(index[1].upper())
    elif index[0].startswith("sciMET_i5"):
        i5_dict[index[0]] = rev_complement(index[1][1:].upper())
    elif index[0].startswith("sciMET_i7"):
        i7_dict[index[0]] = rev_complement(index[1].upper())
# print(tn5_dict)
# print(i5_dict)

## testing for duplicate sequences using set -> should not contain any duplicates
# seq = []
# for key, value in i7_dict.items():
#     seq.append(value)
# unique_seq = set(seq)
# print(len(unique_seq))

with open(output_location + "/barcodes_scimet.txt", "w") as wf:
    i = 0
    for name, seq in i7_dict.items():
        wf.write(name + "\t" + "1" + "\t" + seq + "\n")
        i += 1
        if i >= 9: ## we only used i1-i9
            break


    for name, seq in tn5_dict.items():
        wf.write(name + "\t" + "2" + "\t" + seq + "\n")

    for name, seq in i5_dict.items():
        wf.write(name + "\t" + "3" + "\t" + seq + "\n")
