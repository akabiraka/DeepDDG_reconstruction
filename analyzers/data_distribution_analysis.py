import sys
sys.path.append("../DeepDDG_reconstruction")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Bio.PDB.Polypeptide as Polypeptide

def plot_col_histogram(train_set_df, test_set_df, col_name, 
                       xlabel, ylabel, output_filename, 
                       bins=15, rotate_xticks=False, remove_xticks=False, sort_by_residue=False, save=False):
    train_set_df[col_name].hist(bins=bins, grid=False, label="Train set", color="lightgreen")
    plt.legend()
    test_set_df[col_name].hist(bins=bins, grid=False, label="Test set", color="salmon")
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # locs, labels = plt.xticks()
    # print(locs, labels)
    if rotate_xticks: plt.xticks(rotation=45)
    if remove_xticks: plt.xticks([])
    if sort_by_residue: plt.xticks(np.arange(20), [Polypeptide.index_to_three(i) for i in range(20)], rotation=45)
    if save:
        plt.savefig("output_images/data_distribution/{}.pdf".format(output_filename), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
        plt.close()
    else: plt.show()
    
def parse_mutation_sites(dataset_df, sheet_name, writer_object=None):
    dataset_df["wild_residue"] = dataset_df["mutation"].apply(lambda row: Polypeptide.one_to_three(row.split("_")[0]))
    dataset_df["mutation_site"] = dataset_df["mutation"].apply(lambda row: int(row.split("_")[1]))
    dataset_df["mutant_residue"] = dataset_df["mutation"].apply(lambda row: Polypeptide.one_to_three(row.split("_")[2]))
    if writer_object is None: dataset_df.to_excel("data/dataset_2.xlsx", sheet_name=sheet_name, index=False)
    else: dataset_df.to_excel(writer_object, sheet_name=sheet_name, index=False)
    
    
# parse_mutation_sites(pd.read_excel("data/dataset_1.xlsx", sheet_name="train"), sheet_name="train")
# with pd.ExcelWriter('data/dataset_2.xlsx', mode="a") as writer: 
#     parse_mutation_sites(pd.read_excel("data/dataset_1.xlsx", sheet_name="test"), sheet_name="test", writer_object=writer)
        
train_set_df = pd.read_excel("data/dataset_2.xlsx", sheet_name="train")
test_set_df = pd.read_excel("data/dataset_2.xlsx", sheet_name="test")

# unique_proteins_train = train_set_df["protein_name"].drop_duplicates().to_list() #211s
unique_proteins_train = train_set_df["pdb_id"].drop_duplicates().to_list()
unique_proteins_test = test_set_df["pdb_id"].drop_duplicates().to_list()
print(len(unique_proteins_train), len(unique_proteins_test)) # 209, 37

common_proteins = list(set(unique_proteins_train) & set(unique_proteins_test))
print(len(common_proteins)) # 0

# plot_col_histogram(train_set_df=train_set_df, test_set_df=test_set_df, col_name="pH", 
#                    xlabel="pH", ylabel="Number of mutations", output_filename="pH_vs_mutations_histogram")

# plot_col_histogram(train_set_df=train_set_df, test_set_df=test_set_df, col_name="T", 
#                    xlabel="Temperature (c)", ylabel="Number of mutations", output_filename="temperature_vs_mutations_histogram")

# plot_col_histogram(train_set_df=train_set_df, test_set_df=test_set_df, col_name="ddG", bins=80,
#                    xlabel="ddG (kcal/mol)", ylabel="Number of mutations", output_filename="ddG_vs_mutations_histogram")

# plot_col_histogram(train_set_df=train_set_df, test_set_df=test_set_df, col_name="protein_name", 
#                    xlabel="Proteins", ylabel="Number of mutations", output_filename="proteins_vs_mutations_histogram", 
#                    remove_xticks=True)

# plot_col_histogram(train_set_df=train_set_df, test_set_df=test_set_df, col_name="mutation_site", 
#                    xlabel="Mutation site", ylabel="Number of mutations", output_filename="mutation_site_vs_mutations_histogram")

# plot_col_histogram(train_set_df=train_set_df, test_set_df=test_set_df, col_name="wild_residue", 
#                    xlabel="Amino acids", ylabel="Number of mutations", output_filename="amino_acid_vs_mutations_histogram", 
#                    sort_by_residue=True)

# train_set_df["wild_mutant_residue_pair"] = train_set_df["wild_residue"] + "-" + train_set_df["mutant_residue"]
# test_set_df["wild_mutant_residue_pair"] = test_set_df["wild_residue"] + "-" + test_set_df["mutant_residue"]
# plot_col_histogram(train_set_df=train_set_df, test_set_df=test_set_df, col_name="wild_mutant_residue_pair", 
#                    xlabel="Wild-mutant residue pairs", ylabel="Number of mutations", output_filename="wild_mutant_pair_vs_mutations_histogram",
#                    remove_xticks=True, bins=50)
