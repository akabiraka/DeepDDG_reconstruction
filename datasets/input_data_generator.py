import sys
sys.path.append("../DeepDDG_reconstruction")

import pandas as pd

from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect

# configurations
pdb_dir = "data/pdbs/"
clean_pdb_dir = "data/pdbs_clean/"
CIF = "mmCif"
input_file_path = "data/dataset_2.xlsx"
sheet_name = "train"
n_rows_to_skip = 0
n_rows_to_evalutate = 22

# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)

# data generation
dfs = pd.read_excel(input_file_path, sheet_name=sheet_name)
# print(dfs.columns)

for i, row in dfs.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    pdb_id = row["pdb_id"].lower()
    chain_id = "A"
    PDBData.download_structure(pdb_id=pdb_id)
    clean_wild_protein_structure = PDBData.clean(pdb_id=pdb_id, chain_id=chain_id, selector=ChainAndAminoAcidSelect(chain_id))
    
    print("Protein no: {} -> {}".format(i+1, pdb_id))
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break