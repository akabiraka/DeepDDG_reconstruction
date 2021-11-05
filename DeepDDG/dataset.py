import sys
sys.path.append("../DeepDDG_reconstruction")

import pandas as pd
import numpy as np

import torch
from torch.utils.data import Dataset

class DeepDDGDataset(Dataset):
    def __init__(self, file, device):
        """[summary]

        Args:
            file (str): a dataframe filepath 
        """
        self.data_dir = "data/features_keep/"
        self.mutation_df = pd.read_csv(file)
        self.device = device
        # print(self.mutation_df.shape)
    
    def _parse_a_row(self, row):
        return row["pdb_id"], row["chain_id"], row["mutation"], row["ddG"], row["wild_residue"], row["mutation_site"], row["mutant_residue"]
        
    def __len__(self):
        return self.mutation_df.shape[0]
    
    def __getitem__(self, i):
        pdb_id, chain_id, mutation, ddG, wild_residue, mutation_site, mutant_residue = self._parse_a_row(self.mutation_df.loc[i])
        file_name = pdb_id+"_"+chain_id+"_"+mutation
        
        target_residue_tensor = torch.load(self.data_dir+file_name+".pt").unsqueeze(dim=0)
        all_neighbor_residue_tensor = torch.load(self.data_dir+file_name+"_neighbors.pt")
        ddG = torch.tensor(ddG, dtype=torch.float32).unsqueeze(dim=0)
        
        broadcast_shape = (all_neighbor_residue_tensor.shape[0], target_residue_tensor.shape[1])
        target_residue_tensor = torch.broadcast_to(target_residue_tensor, broadcast_shape)
        pair_tensors = torch.cat((target_residue_tensor, all_neighbor_residue_tensor), dim=1)
        # print(target_residue_tensor.dtype, target_residue_tensor.shape)
        # print(all_neighbor_residue_tensor.dtype, all_neighbor_residue_tensor.shape)
        # print(ddG.dtype, ddG.shape, ddG)
        # print(pair_tensors.shape, pair_tensors.device)
        
        return pair_tensors, ddG
    
    
# train_dataset = DeepDDGDataset(file="data/dataset_4.xlsx", device="cuda")
# print(train_dataset.__len__())
# ddG = train_dataset.__getitem__(i=0)
# target_residue_tensor, all_neighbor_residue_tensor, ddG = train_dataset.__getitem__(i=0)
