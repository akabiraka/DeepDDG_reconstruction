import sys
sys.path.append("../DeepDDG_reconstruction")
from Bio.PDB import PDBParser
def get_starting_residue_index(pdb_file):
    with open(pdb_file) as pdb_reader:
        for line in pdb_reader:
            if line.split()[5].isnumeric(): 
                return int(line.split()[5])

# clean_pdb_file = "data/pdbs_clean/1lmb1.pdb" 
# print(get_starting_residue_index(pdb_file=clean_pdb_file))

# for model in PDBParser(QUIET=True).get_structure("", clean_pdb_file):
#     print(model)