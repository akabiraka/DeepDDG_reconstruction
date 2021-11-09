import pandas as pd
from Bio.PDB import Polypeptide

def cleate_mutation_columns(dataset_df):
    dataset_df["wild_residue"] = dataset_df["mutation"].apply(lambda row: Polypeptide.one_to_three(row.split("_")[0]))
    dataset_df["mutation_site"] = dataset_df["mutation"].apply(lambda row: int(row.split("_")[1]))
    dataset_df["mutant_residue"] = dataset_df["mutation"].apply(lambda row: Polypeptide.one_to_three(row.split("_")[2]))
    return dataset_df

train_set_df = pd.read_excel("data/dataset_1.xlsx", sheet_name="train")
test_set_df = pd.read_excel("data/dataset_1.xlsx", sheet_name="test")

# removing duplicate entries for train set by grouping onto pdb_id and mutation
train_set_df = train_set_df[["pdb_id", "mutation", "ddG"]] 
grouped_train_set_df = train_set_df.groupby(by=["pdb_id", "mutation"], as_index=False).mean()
# print(grouped_train_set_df)
test_set_df = test_set_df[["pdb_id", "mutation", "ddG"]]
grouped_test_set_df = test_set_df.groupby(by=["pdb_id", "mutation"], as_index=False).mean()
# print(grouped_test_set_df)


unique_proteins_train = grouped_train_set_df["pdb_id"].drop_duplicates().to_list()
unique_proteins_test = grouped_test_set_df["pdb_id"].drop_duplicates().to_list()
print(len(unique_proteins_train), len(unique_proteins_test)) # 209, 37
common_proteins = list(set(unique_proteins_train) & set(unique_proteins_test))
print(len(common_proteins)) # 0

# separating mutation into 3 columns
grouped_train_set_df = cleate_mutation_columns(grouped_train_set_df)
grouped_test_set_df = cleate_mutation_columns(grouped_test_set_df)

# saving data
grouped_train_set_df.to_excel("data/dataset_3_train.xlsx", index=False)
grouped_test_set_df.to_excel('data/dataset_3_test.xlsx', index=False)