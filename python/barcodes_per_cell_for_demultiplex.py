## this script takes the csv file containing the barcode sequences and make an excel file containing cell types for assignment
## arguments needed: csv file containing the index name and sequence; output location
## the sequences should be the reverse complement of the sequences in the csv file

import sys
import pandas as pd
import itertools

## import csv file
csv_file = sys.argv[1]
output_location = sys.argv[2]

pd.set_option("display.max_rows", 400) # this sets up the number of rows displayed
pd.set_option("display.max_columns", 10)
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

i7_dict_scimet_n9 = dict(itertools.islice(i7_dict.items(), 0, 3))
i7_dict_scimet_mg = dict(itertools.islice(i7_dict.items(), 3, 6))
i7_dict_scimet_em = dict(itertools.islice(i7_dict.items(), 6, 9))

# print(i7_dict_scimet_n9)
# print(i7_dict_scimet_mg)
# print(i7_dict_scimet_em)

## combination of sequences from i7, Tn5, i5
combinations_n9 = {}
combinations_mg = {}
combinations_em = {}

for i7_keys, i7_values in i7_dict_scimet_n9.items():
    for met_keys, met_values in tn5_dict.items():
        for i5_keys, i5_values in i5_dict.items():
            combinations_n9[i7_keys + "_" + met_keys + "_" + i5_keys] = i7_values + met_values + i5_values

for i7_keys, i7_values in i7_dict_scimet_mg.items():
    for met_keys, met_values in tn5_dict.items():
        for i5_keys, i5_values in i5_dict.items():
            combinations_mg[i7_keys + "_" + met_keys + "_" + i5_keys] = i7_values + met_values + i5_values

for i7_keys, i7_values in i7_dict_scimet_em.items():
    for met_keys, met_values in tn5_dict.items():
        for i5_keys, i5_values in i5_dict.items():
            combinations_em[i7_keys + "_" + met_keys + "_" + i5_keys] = i7_values + met_values + i5_values

## create dataframe from dict + add "Experiment_name" column
df_n9 = pd.DataFrame(combinations_n9.items(), columns = ["seqName", "Sequence"])
df_mg = pd.DataFrame(combinations_mg.items(), columns = ["seqName", "Sequence"])
df_em = pd.DataFrame(combinations_em.items(), columns = ["seqName", "Sequence"])

# df_n9 = df_n9.insert(0, "Experiment_name", range(0, len(df_n9))
df_n9["Cell_num"] = range(1, len(df_n9)+1); df_n9["Exp_name"] = "sciMETN9"; df_n9["Experiment_name"] = df_n9["Exp_name"] + "_" + df_n9["Cell_num"].astype(str); exp_name = df_n9["Experiment_name"]; df_n9 = df_n9.drop(["Cell_num", "Exp_name", "Experiment_name"], axis = 1); df_n9.insert(0, "Experiment_name", exp_name)
df_mg["Cell_num"] = range(1, len(df_mg)+1); df_mg["Exp_name"] = "sciMETMG"; df_mg["Experiment_name"] = df_mg["Exp_name"] + "_" + df_mg["Cell_num"].astype(str); exp_name = df_mg["Experiment_name"]; df_mg = df_mg.drop(["Cell_num", "Exp_name", "Experiment_name"], axis = 1); df_mg.insert(0, "Experiment_name", exp_name)
df_em["Cell_num"] = range(1, len(df_em)+1); df_em["Exp_name"] = "sciMETEM"; df_em["Experiment_name"] = df_em["Exp_name"] + "_" + df_em["Cell_num"].astype(str); exp_name = df_em["Experiment_name"]; df_em = df_em.drop(["Cell_num", "Exp_name", "Experiment_name"], axis = 1); df_em.insert(0, "Experiment_name", exp_name)

## save the dataframe to excel file
with pd.ExcelWriter('single_cell_types_for_demultiplex.xlsx') as writer:
    df_n9.to_excel(writer, sheet_name='df_n9')
    df_mg.to_excel(writer, sheet_name='df_mg')
    df_em.to_excel(writer, sheet_name='df_em')
