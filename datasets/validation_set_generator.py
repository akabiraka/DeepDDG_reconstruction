import pandas as pd
import random

def generate_validation_set(n_pdbids=20):
    train_set_df = pd.read_csv("data/dataset_4_train_keep.csv")

    unique_train_pdbids = train_set_df["pdb_id"].drop_duplicates().to_list()
    unique_val_pdbids = random.sample(unique_train_pdbids, n_pdbids)

    columns = ["pdb_id", "chain_id", "mutation", "ddG", "wild_residue", "mutation_site", "mutant_residue"]
    val_set_df = pd.DataFrame(columns=columns)
    for val_pdbid in unique_val_pdbids:
        val_data = train_set_df[train_set_df["pdb_id"]==val_pdbid]
        train_set_df.drop(val_data.index, inplace=True)
        val_set_df = val_set_df.append(val_data)

    print(train_set_df.shape, val_set_df.shape)
    # print(val_set_df)
    val_set_df.to_csv("data/dataset_5_validation.csv", index=False)
    train_set_df.to_csv("data/dataset_5_train.csv", index=False)