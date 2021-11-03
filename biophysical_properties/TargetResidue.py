import re
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

class TargetResidue(object):
    def __init__(self) -> None:
        super().__init__()
        self.backbone_dihedral = BackboneDihedral()
        self.sasa = SASA()
        self.secondary_structure = SecondaryStructure() 
        self.pssm = PSSM()
        self.hydrogen_bond = HydrogenBond()
        self.distance_and_orientation = DistanceAndOrientation()

    def get_features(self, clean_pdb_file, wild_fasta_file, mutant_fasta_file, wild_residue, 
                     mutant_residue, chain_id, mutation_site, starting_residue_id):
        angles = self.backbone_dihedral.of_a_residue(pdb_file=clean_pdb_file, residue_num=mutation_site, chain_id=chain_id, return_type="both")

        self.sasa.set_up(pdb_file=clean_pdb_file)
        sasa_value = self.sasa.of_a_residue(residue_index=mutation_site)

        ss_one_hot = self.secondary_structure.of_a_residue(pdb_file=clean_pdb_file, residue_index=mutation_site-starting_residue_id, return_type="one-hot")

        self.pssm.set_up(wild_fasta_file)
        wild_residue_pssm_score = self.pssm.of_a_residue(residue_index=mutation_site-starting_residue_id)
        self.pssm.set_up(mutant_fasta_file)
        mutant_residue_pssm_score = self.pssm.of_a_residue(residue_index=mutation_site-starting_residue_id)

        wild_residue_type = np.array(list('{0:05b}'.format(Polypeptide.three_to_index(wild_residue))), dtype=np.float32)
        mutant_residue_type = np.array(list('{0:05b}'.format(Polypeptide.three_to_index(mutant_residue))), dtype=np.float32)
        
        # concatenating all target residue features
        target_residue_features = np.concatenate((angles, sasa_value, ss_one_hot, wild_residue_pssm_score, mutant_residue_pssm_score, wild_residue_type, mutant_residue_type))
        # print(target_residue_features.shape, target_residue_features.dtype, target_residue_features)
        # print(angles, sasa_value, ss_one_hot, wild_residue_pssm_score, mutant_residue_pssm_score, wild_residue_type, mutant_residue_type)
        return target_residue_features


# clean_pdb_file = "data/pdbs_clean/1a43A.pdb"
# wild_fasta_file = "data/fastas/1a43A.fasta"  
# mutant_fasta_file = "data/fastas/1a43A_mutant.fasta"
# wild_residue = "TRP"
# mutant_residue = "ALA"
# chain_id = "A"
# mutation_site = 184
# starting_residue_id = pdb_utils.get_starting_residue_index(pdb_file=clean_pdb_file)

# TR = TargetResidue()
# TR.get_features(clean_pdb_file=clean_pdb_file, wild_fasta_file=wild_fasta_file, 
#                 mutant_fasta_file=mutant_fasta_file, wild_residue=wild_residue, 
#                 mutant_residue=mutant_residue, chain_id=chain_id, 
#                 mutation_site=mutation_site, starting_residue_id=starting_residue_id)