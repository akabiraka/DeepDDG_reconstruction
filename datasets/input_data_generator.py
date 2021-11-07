from os import write
import sys
sys.path.append("../DeepDDG_reconstruction")

import csv
import pandas as pd
import numpy as np
import torch
from torch import nn
import torch.optim as optim
from Bio.PDB import Polypeptide
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
# input_file_path = "data/bad_things_check.xlsx"
input_file_path = "data/dataset_3_train.xlsx"
output_file_path = "data/dataset_4_train.csv"
n_rows_to_skip = 203
n_rows_to_evalutate = 10000
N_neighbors = 15

# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)
target_residue = TargetResidue()
neighbor_residue = NeighborResidue()
cleanslate = CleanSlate()
cleanslate.clean_all()


# data generation
dfs = pd.read_excel(input_file_path)
# print(dfs.columns)

def validate_chain_id(given_pdb_id):
    pdb_id = row["pdb_id"].lower()[:4]
    
    if pdb_id == "1lmb": return "4"
    if pdb_id == "1tup": return "A"
    
    chain_id = PDBData.get_first_chain_id(pdb_id=pdb_id)
    if len(given_pdb_id) > 7 and given_pdb_id[4]==":" and given_pdb_id[6]==":":
        chain_id = given_pdb_id[5].upper()
    return chain_id

columns = ["pdb_id", "chain_id", "mutation", "ddG", "wild_residue", "mutation_site", "mutant_residue"]
with open(output_file_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        
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
    row=[pdb_id, chain_id, mutation, row["ddG"], wild_residue, mutation_site, mutant_residue]
    with open(output_file_path, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(row)
    
    clean_pdb_file = pdbs_clean_dir+pdb_id+chain_id+".pdb"
    clean_wild_protein_structure = PDBData.clean(pdb_id=pdb_id, chain_id=chain_id, selector=ChainAndAminoAcidSelect(chain_id))
    wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
    mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_"+mutation+".fasta"
    PDBData.generate_fasta_from_pdb(pdb_id, chain_id, clean_pdb_file, save_as_fasta=True, output_fasta_dir="data/fastas/")
    PDBData.create_mutant_fasta_file(wild_fasta_file, mutant_fasta_file, mutation_site, wild_residue)
    residue_ids_dict = PDBData.get_residue_ids_dict(pdb_file=clean_pdb_file, chain_id="A")
    zero_based_mutation_site = residue_ids_dict.get(mutation_site)
    print("Row no:{}->{}{}, mutation:{}, zero_based_mutation_site:{}".format(i+1, pdb_id, chain_id, mutation, zero_based_mutation_site))
    
    
    # computing target residue features
    target_residue_features = target_residue.get_features(clean_pdb_file=clean_pdb_file, wild_fasta_file=wild_fasta_file, 
                                mutant_fasta_file=mutant_fasta_file, wild_residue=wild_residue, 
                                mutant_residue=mutant_residue, chain_id=chain_id, 
                                mutation_site=mutation_site, zero_based_mutation_site=zero_based_mutation_site)
    # print(target_residue_features.shape, target_residue_features.dtype, target_residue_features)
    
    # computing neighbor residue features
    all_neighbor_features = []
    n_neighbor_residue_ids = neighbor_residue.get_n_neighbor_residue_ids(pdb_file=clean_pdb_file, 
                                                                         chain_id=chain_id, 
                                                                         center_residue_id=mutation_site, N=N_neighbors)
    print("All neighbors: ", n_neighbor_residue_ids)
    for neighbor_residue_id in n_neighbor_residue_ids:
        print("Neighbor residue id: ", neighbor_residue_id)
        neighbor_residue_features = neighbor_residue.get_features(clean_pdb_file=clean_pdb_file, chain_id=chain_id, 
                                      mutation_site=mutation_site, zero_based_mutation_site=zero_based_mutation_site, 
                                      neighbor_residue_id=neighbor_residue_id)
        # print(neighbor_residue_features.shape, neighbor_residue_features.dtype, neighbor_residue_features)
        all_neighbor_features.append(neighbor_residue_features)
        
    # print(np.array(all_neighbor_features).shape)    
    
    file_name = pdb_id+"_"+chain_id+"_"+mutation
    torch.save(torch.tensor(np.array(target_residue_features, dtype=np.float32)), "data/features/"+file_name+".pt")
    torch.save(torch.tensor(np.array(all_neighbor_features, dtype=np.float32)), "data/features/"+file_name+"_neighbors.pt")
    
    target_residue_tensor = torch.load("data/features/"+file_name+".pt")
    print(target_residue_tensor.size())
    all_neighbor_residue_tensor = torch.load("data/features/"+file_name+"_neighbors.pt")
    print(all_neighbor_residue_tensor.size())

    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break
    
# updated_dataset.to_excel("data/dataset_4.xlsx", index=False)