import pandas as pd

import torch
from torch.utils.data import Dataset

class DeepDDGDataset(Dataset):
    def __init__(self, file):
        """[summary]

        Args:
            file (str): a dataframe filepath 
        """
        self.mutation_df = pd.read_excel(file)
        # print(self.mutation_df.shape)
    
    def _parse_a_row(self, row):
        return row["pdb_id"], row["chain_id"], row["mutation"], row["ddG"], row["wild_residue"], row["mutation_site"], row["mutant_residue"]
        
    def __len__(self):
        return self.mutation_df.shape[0]
    
    def __getitem__(self, i):
        pdb_id, chain_id, mutation, ddG, wild_residue, mutation_site, mutant_residue = self._parse_a_row(self.mutation_df.loc[i])
        file_name = pdb_id+"_"+chain_id+"_"+mutation
        
        target_residue_tensor = torch.load("data/features/"+file_name+".pt")
        all_neighbor_residue_tensor = torch.load("data/features/"+file_name+"_neighbors.pt")
        
        return target_residue_tensor, all_neighbor_residue_tensor
    
    
# train_dataset = DeepDDGDataset(file="data/dataset_4.xlsx")
# print(train_dataset.__len__())
# target_residue_tensor, all_neighbor_residue_tensor = train_dataset.__getitem__(i=0)
# print(target_residue_tensor.dtype, target_residue_tensor.shape)
# print(all_neighbor_residue_tensor.dtype, all_neighbor_residue_tensor.shape)