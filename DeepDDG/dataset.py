import sys

from torch.utils import data
sys.path.append("../DeepDDG_reconstruction")

import pandas as pd
import numpy as np

import torch
from torch.utils.data import Dataset

class DeepDDGDataset(Dataset):
    def __init__(self, file, data_dir=None, do_normalize=False):
        """[summary]

        Args:
            file (str): a dataframe filepath 
        """
        self.do_normalize = do_normalize
        self.data_dir = "data/features/" if data_dir==None else data_dir
        self.mutation_df = pd.read_csv(file)
        self.mutation_df = self.normalize(self.mutation_df, -1, 1)
        # print(self.mutation_df.shape)
    
    def _parse_a_row(self, row):
        class_label = 1.0 if row["ddG"]>=0 else 0.0
        return row["pdb_id"], row["chain_id"], row["mutation"], row["ddG"], row["wild_residue"], \
            row["mutation_site"], row["mutant_residue"], class_label, row["normalized_ddg"]
        
    def __len__(self):
        return self.mutation_df.shape[0]
    
    def __getitem__(self, i):
        pdb_id, chain_id, mutation, ddG, wild_residue, mutation_site, mutant_residue, class_label, normalized_ddg = self._parse_a_row(self.mutation_df.loc[i])
        file_name = pdb_id+"_"+chain_id+"_"+mutation
        
        target_residue_tensor = torch.load(self.data_dir+file_name+".pt").unsqueeze(dim=0)
        all_neighbor_residue_tensor = torch.load(self.data_dir+file_name+"_neighbors.pt")
        all_neighbor_residue_tensor[:, 14] = all_neighbor_residue_tensor[:, 14]/20 # dividing Ca-Ca distance value by 20A, suggested by the paper
        normalized_ddg = torch.tensor(normalized_ddg, dtype=torch.float32).unsqueeze(dim=0)
        ddG = torch.tensor(ddG, dtype=torch.float32).unsqueeze(dim=0)
        class_label = torch.tensor(class_label, dtype=torch.float32).unsqueeze(dim=0)
        
        broadcast_shape = (all_neighbor_residue_tensor.shape[0], target_residue_tensor.shape[1])
        target_residue_tensor = torch.broadcast_to(target_residue_tensor, broadcast_shape)
        pair_tensors = torch.cat((target_residue_tensor, all_neighbor_residue_tensor), dim=1)
        # print(target_residue_tensor.dtype, target_residue_tensor.shape)
        # print(all_neighbor_residue_tensor.dtype, all_neighbor_residue_tensor.shape)
        # print(ddG.dtype, ddG.shape, ddG)
        # print(pair_tensors.shape, pair_tensors.device)
        
        if self.do_normalize:
            return pair_tensors, normalized_ddg, class_label
        return pair_tensors, ddG, class_label
    
    def normalize(self, df, x, y):
        # x, y = -1, 1
        # df = pd.read_csv("data/dataset_5_train.csv")
        all_ddgs = df["ddG"].to_list()
        min_ddg = min(all_ddgs)
        max_ddg = max(all_ddgs)
        range_ddg = max_ddg - min_ddg
        # print(max_ddg, min_ddg, range_ddg)

        for i, row in df.iterrows():
            ddg = row["ddG"]
            ddg = (ddg - min_ddg)/range_ddg
            range_scale = y-x
            normalized_ddg = (ddg*range_scale)+x
            df.loc[i, "normalized_ddg"] = normalized_ddg
        
        return df    
        # df["normalized_ddg"].hist(bins=35, grid=False, label="Train set", color="lightgreen")  
        # plt.show() 
# train_dataset = DeepDDGDataset(file="data/dataset_4_train_keep.csv", data_dir="data/features_train/")
# print(train_dataset.__len__())
# train_dataset.__getitem__(i=0)
# target_residue_tensor, all_neighbor_residue_tensor, ddG = train_dataset.__getitem__(i=0)
