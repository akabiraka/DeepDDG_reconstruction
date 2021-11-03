import sys
sys.path.append("../DeepDDG_reconstruction")

import numpy as np
from Bio.PDB import PDBParser, Polypeptide
from utils import pdb_utils
from biophysical_properties.BackboneDihedral import BackboneDihedral
from biophysical_properties.SASA import SASA
from biophysical_properties.SecondaryStructure import SecondaryStructure
from biophysical_properties.PSSM import PSSM
from biophysical_properties.DistanceAndOrientation import DistanceAndOrientation
from biophysical_properties.HydrogenBond import HydrogenBond

class NeighborResidue(object):
    def __init__(self) -> None:
        super().__init__()
        self.backbone_dihedral = BackboneDihedral()
        self.sasa = SASA()
        self.secondary_structure = SecondaryStructure() 
        self.pssm = PSSM()
        self.hydrogen_bond = HydrogenBond()
        self.distance_and_orientation = DistanceAndOrientation()


    def get_n_neighbor_residue_ids(self, pdb_file, chain_id, center_residue_id, N):
        pdb_id = pdb_file.split("/")[2].split(".")[0]
        center_residue = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)[0][chain_id][center_residue_id]
        residues = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)[0][chain_id].get_residues()
        
        starting_residue_index = pdb_utils.get_starting_residue_index(pdb_file=pdb_file)

        residue_id_vs_distance = []
        for i, residue in enumerate(residues, starting_residue_index):
            if i==center_residue_id: continue
            diff_vector = center_residue["CA"].coord - residue["CA"].coord
            distance = np.sqrt(np.sum(diff_vector * diff_vector))
            residue_id_vs_distance.append([i, distance])
        
        n_neighbor_residue_ids = np.array(sorted(residue_id_vs_distance, key=lambda x: x[1]))[:N, 0]
        return n_neighbor_residue_ids.astype(int).tolist()
            

    def get_features(self, clean_pdb_file, chain_id, mutation_site, starting_residue_id, neighbor_residue_id):
        neighbor_residue = PDBParser(QUIET=True).get_structure("", clean_pdb_file)[0][chain_id][neighbor_residue_id]

        angles = self.backbone_dihedral.of_a_residue(pdb_file=clean_pdb_file, residue_num=neighbor_residue_id, chain_id=chain_id, return_type="both")

        self.sasa.set_up(pdb_file=clean_pdb_file)
        sasa_value = self.sasa.of_a_residue(residue_index=neighbor_residue_id)

        ss_one_hot = self.secondary_structure.of_a_residue(pdb_file=clean_pdb_file, residue_index=mutation_site-starting_residue_id, return_type="one-hot")
        
        self.hydrogen_bond.set_up(pdb_file=clean_pdb_file)
        num_of_hydrogen_bonds = self.hydrogen_bond.get(target_residue_id=mutation_site, neighbor_residue_id=neighbor_residue_id)
   
        ca_ca_distance, ca_ca_unit_vector, ca_c_unit_vector, ca_n_unit_vector = self.distance_and_orientation.get(pdb_file=clean_pdb_file, chain_id=chain_id, target_residue_id=mutation_site, neighbor_residue_id=neighbor_residue_id)
    
        neighbor_residue_type = np.array(list('{0:05b}'.format(Polypeptide.three_to_index(neighbor_residue.get_resname()))), dtype=np.float32)
       

        # concatenating all neighbor residue features
        neighbor_residue_features = np.concatenate((angles, sasa_value, ss_one_hot, num_of_hydrogen_bonds, np.array([ca_ca_distance]), ca_ca_unit_vector, ca_c_unit_vector, ca_n_unit_vector, neighbor_residue_type))
        # print(neighbor_residue_features.shape, neighbor_residue_features.dtype, neighbor_residue_features)
        return neighbor_residue_features




# clean_pdb_file = "data/pdbs_clean/1a43A.pdb" 
# chain_id = "A"
# mutation_site = 184
# starting_residue_id = pdb_utils.get_starting_residue_index(pdb_file=clean_pdb_file)

# NR = NeighborResidue()

# n_neighbor_residue_ids = NR.get_n_neighbor_residue_ids(pdb_file=clean_pdb_file, chain_id=chain_id, center_residue_id=mutation_site, N=5)
# print(n_neighbor_residue_ids)

# for neighbor_residue_id in n_neighbor_residue_ids:
#     NR.get_features(clean_pdb_file, chain_id, mutation_site, starting_residue_id, neighbor_residue_id)