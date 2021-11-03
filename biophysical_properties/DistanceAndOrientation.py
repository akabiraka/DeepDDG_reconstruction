import sys
sys.path.append("../DeepDDG_reconstruction")

import numpy as np
from Bio.PDB import PDBParser

class DistanceAndOrientation(object):
    def __init__(self) -> None:
        super().__init__()

    def get(self, pdb_file, chain_id, target_residue_id, neighbor_residue_id):
        target_residue = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id][target_residue_id]
        neighbor_residue = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id][neighbor_residue_id]
        
        ca_ca_diff_vector = target_residue["CA"].coord - neighbor_residue["CA"].coord
        ca_ca_distance = np.sqrt(np.sum(ca_ca_diff_vector * ca_ca_diff_vector))
        
        ca_ca_unit_vector = ca_ca_diff_vector / np.linalg.norm(ca_ca_diff_vector)

        ca_c_diff_vector = neighbor_residue["CA"].coord - neighbor_residue["C"].coord
        ca_c_unit_vector = ca_c_diff_vector / np.linalg.norm(ca_c_diff_vector)

        ca_n_diff_vector = neighbor_residue["CA"].coord - neighbor_residue["N"].coord
        ca_n_unit_vector = ca_n_diff_vector / np.linalg.norm(ca_n_diff_vector)

        return ca_ca_distance, ca_ca_unit_vector, ca_c_unit_vector, ca_n_unit_vector
    