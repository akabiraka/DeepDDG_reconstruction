import sys
sys.path.append("../DeepDDG_reconstruction")

import csv
import pandas as pd
import os
from utils import pdb_utils
from objects.PDBData import PDBData
from biophysical_properties.PSSM import PSSM
from objects.Selector import ChainAndAminoAcidSelect

def validate_chain_id(given_pdb_id):
    pdb_id = row["pdb_id"].lower()[:4]
    
    if pdb_id == "1lmb": return "4"
    if pdb_id == "1tup": return "A"
    
    chain_id = PDBData.get_first_chain_id(pdb_id=pdb_id)
    if len(given_pdb_id) > 7 and given_pdb_id[4]==":" and given_pdb_id[6]==":":
        chain_id = given_pdb_id[5].upper()
    return chain_id



# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"
input_file_path = "data/dataset_3_train.xlsx"
# input_file_path = "data/bad_things_check.xlsx"
output_file_path = "data/dataset_4_train.csv"
n_rows_to_skip = 0
n_rows_to_evalutate = 2
N_neighbors = 15

columns = ["pdb_id", "chain_id", "mutation", "ddG", "wild_residue", "mutation_site", "mutant_residue"]
with open(output_file_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        
# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)
pssm = PSSM()

# data generation
dfs = pd.read_excel(input_file_path)

    
# i = int(os.environ["SLURM_ARRAY_TASK_ID"])    
i=208
unique_pdb_ids = dfs["pdb_id"].drop_duplicates().to_list()
unique_pdb_ids.sort()
ith_pdb_id = unique_pdb_ids[i]
ith_protein_mutation_dfs = dfs[dfs["pdb_id"]==ith_pdb_id]
# print(ith_protein_mutation_dfs)

for i, row in ith_protein_mutation_dfs.iterrows():
    # if i+1 <= n_rows_to_skip: return
    
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
    starting_residue_id = pdb_utils.get_first_residue_id(pdb_file=clean_pdb_file, chain_id=chain_id)
    zero_based_mutation_site = mutation_site-starting_residue_id
    print("Row no:{}->{}{}, mutation:{}, first_residue_id:{}, zero_based_mutation_site:{}".format(i+1, pdb_id, chain_id, mutation, starting_residue_id, zero_based_mutation_site))


    pssm.set_up(wild_fasta_file)
    pssm.set_up(mutant_fasta_file)