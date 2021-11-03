import sys
sys.path.append("../DeepDDG_reconstruction")
from Bio.PDB import PDBParser
# def get_starting_residue_index(pdb_file):
#     with open(pdb_file) as pdb_reader:
#         for line in pdb_reader:
#             if line.split()[5].isnumeric(): 
#                 return int(line.split()[5])

def get_first_residue_id(pdb_file, chain_id):
    residue_list = list(PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id].get_residues())
    _, residue_id, _ = residue_list[0].id
    return residue_id
            
def get_last_residue_id(pdb_file, chain_id):
    residue_list = list(PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id].get_residues())
    _, residue_id, _ = residue_list[len(residue_list)-1].id
    return residue_id



# clean_pdb_file = "data/pdbs_clean/1lmb1.pdb" 
# print(get_starting_residue_index(pdb_file=clean_pdb_file))

# print(get_last_residue_id("data/pdbs_clean/1a43A.pdb", "A"))
# print(get_first_residue_id("data/pdbs_clean/1a43A.pdb", "A"))

# PDBParser(QUIET=True).get_structure("", "data/pdbs_clean/1a7cA.pdb")[0]["A"][336]