import sys
sys.path.append("../DeepDDG_reconstruction")

import pandas as pd
import numpy as np
import torch
from torch import nn
import torch.optim as optim
from Bio.PDB import Polypeptide
from utils import pdb_utils
from utils.CleanSlate import CleanSlate
from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
from biophysical_properties.TargetResidue import TargetResidue
from biophysical_properties.NeighborResidue import NeighborResidue

# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"
input_file_path = "data/dataset_2.xlsx"
sheet_name = "train"
n_rows_to_skip = 5444
n_rows_to_evalutate = 2000
N_neighbors = 5

# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)
target_residue = TargetResidue()
neighbor_residue = NeighborResidue()
cleanslate = CleanSlate()
cleanslate.clean_all()


# data generation
dfs = pd.read_excel(input_file_path, sheet_name=sheet_name)
# print(dfs.columns)

def validate_chain_id(given_pdb_id):
    pdb_id = row["pdb_id"].lower()[:4]
    
    if pdb_id == "1lmb": return "4"
    if pdb_id == "1tup": return "A"
    
    chain_id = PDBData.get_first_chain_id(pdb_id=pdb_id)
    if len(given_pdb_id) > 7 and given_pdb_id[4]==":" and given_pdb_id[6]==":":
        chain_id = given_pdb_id[5].upper()
    return chain_id

for i, row in dfs.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    pdb_id = row["pdb_id"].lower()[:4]
    if pdb_id=="2a01": continue
    
    PDBData.download_structure(pdb_id=pdb_id)
    
    chain_id = validate_chain_id(row["pdb_id"])
    mutation = row["mutation"]
    mutation_site = int(row["mutation_site"])
    wild_residue = row["wild_residue"]
    mutant_residue = row["mutant_residue"]
    # print(pdb_id, chain_id)
    
    clean_pdb_file = pdbs_clean_dir+pdb_id+chain_id+".pdb"
    clean_wild_protein_structure = PDBData.clean(pdb_id=pdb_id, chain_id=chain_id, selector=ChainAndAminoAcidSelect(chain_id))
    wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
    mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_mutant.fasta"
    PDBData.generate_fasta_from_pdb(pdb_id, chain_id, clean_pdb_file, save_as_fasta=True, output_fasta_dir="data/fastas/")
    PDBData.create_mutant_fasta_file(wild_fasta_file, fastas_dir, mutation, mutation_site, wild_residue)
    starting_residue_id = pdb_utils.get_starting_residue_index(pdb_file=clean_pdb_file)
    zero_based_mutation_site = mutation_site-starting_residue_id
    print("Protein no: {} -> {}, {}, {}, {}, {}".format(i+1, pdb_id, chain_id, mutation, starting_residue_id, zero_based_mutation_site))
    
    
    # # computing target residue features
    # target_residue_features = target_residue.get_features(clean_pdb_file=clean_pdb_file, wild_fasta_file=wild_fasta_file, 
    #                             mutant_fasta_file=mutant_fasta_file, wild_residue=wild_residue, 
    #                             mutant_residue=mutant_residue, chain_id=chain_id, 
    #                             mutation_site=mutation_site, starting_residue_id=starting_residue_id)
    # # print(target_residue_features.shape, target_residue_features.dtype, target_residue_features)
    
    # # computing neighbor residue features
    # all_neighbor_features = []
    # n_neighbor_residue_ids = neighbor_residue.get_n_neighbor_residue_ids(pdb_file=clean_pdb_file, 
    #                                                                      chain_id=chain_id, 
    #                                                                      center_residue_id=mutation_site, N=N_neighbors)
    # for neighbor_residue_id in n_neighbor_residue_ids:
    #     neighbor_residue_features = neighbor_residue.get_features(clean_pdb_file=clean_pdb_file, chain_id=chain_id, 
    #                                   mutation_site=mutation_site, starting_residue_id=starting_residue_id, 
    #                                   neighbor_residue_id=neighbor_residue_id)
    #     # print(neighbor_residue_features.shape, neighbor_residue_features.dtype, neighbor_residue_features)
    #     all_neighbor_features.append(neighbor_residue_features)
        
    # print(np.array(all_neighbor_features).shape)    
    
    # file_name = pdb_id+"_"+chain_id+"_"+wild_residue+"_"+str(mutation_site)+"_"+mutant_residue
    # torch.save(torch.tensor(np.array(target_residue_features)), "data/features/"+file_name+".pt")
    # torch.save(torch.tensor(np.array(all_neighbor_features)), "data/features/"+file_name+"_neighbors.pt")
    
    # target_residue_tensor = torch.load("data/features/"+file_name+".pt")
    # all_neighbor_residue_tensor = torch.load("data/features/"+file_name+"_neighbors.pt")
    # print(target_residue_tensor.size(), all_neighbor_residue_tensor.size())
    
    # # model compatibility test with data
    # from models.DeepDDG import *
    
    # device = 'cuda' if torch.cuda.is_available() else 'cpu'
    # print('Using {} device'.format(device))
    
    # srp_model = SRP(in_features=51, out_features=20).to(device)
    # deepddg_model = DeepDDG(in_features=N_neighbors*20).to(device)
    
    # criterion = nn.MSELoss()
    # srp_optimizer = optim.Adam(srp_model.parameters(), lr=0.0008)
    # deepddg_optimizer = optim.Adam(deepddg_model.parameters(), lr=0.0008)
    
    # # running the model
    # srp_outs = []
    # for neighbor_residue_tensor in all_neighbor_residue_tensor:
    #     pair_tensor = torch.cat((target_residue_tensor, neighbor_residue_tensor))
    #     pair_tensor = pair_tensor.to(device=device, dtype=torch.float32)
    #     pair_tensor.unsqueeze_(dim=0)
    #     print(pair_tensor.shape, pair_tensor.device)
    #     srp_outs.append(srp_model(pair_tensor))
        
    # concatenated_srp_outs = torch.cat(srp_outs, dim=1)
    # print(concatenated_srp_outs.shape)
    # ddg_pred = deepddg_model(concatenated_srp_outs)
    # print(ddg_pred.shape, ddg_pred)
    
    # # computing loss, backpropagate and optimizing model
    # ddg_target = torch.randn(1, 1, device=device)
    # loss = criterion(ddg_target, ddg_pred)
    # print(loss)
    # loss.backward()
    # srp_optimizer.step()
    # deepddg_optimizer.step()
        
    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break