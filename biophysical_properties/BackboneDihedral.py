import sys
sys.path.append("../DeepDDG_reconstruction")

import math
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral

class BackboneDihedral(object):
    def __init__(self) -> None:
        super().__init__()
        
    def get_angles(self, pdb_file, residue_num, chain_id="A"):
        """
        # used the following two links to compute phi, psi and omega angels
            # https://biopython.org/docs/dev/api/Bio.PDB.internal_coords.html#Bio.PDB.internal_coords.IC_Residue.pick_angle
            # https://biopython.org/docs/latest/api/Bio.PDB.vectors.html?highlight=calc_dihedral#Bio.PDB.vectors.calc_dihedral
            # (sN, sCA, sC, nN)   # psi
            # (sCA, sC, nN, nCA)  # omega i+1
            # (sC, nN, nCA, nC)   # phi i+1

        Args:
            pdb_file ([type]): [description]
            residue_num ([type]): [description]
            chain_id (str, optional): [description]. Defaults to "A".
        Returns:
            phi, psi and omega 
        """
        residue = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id][residue_num]
        next_residue = PDBParser(QUIET=True).get_structure("", pdb_file)[0]["A"][residue_num+1]
        
        sN, sCA, sC = residue["N"].get_vector(), residue["CA"].get_vector(), residue["C"].get_vector()
        nN, nCA, nC = next_residue["N"].get_vector(), next_residue["CA"].get_vector(), next_residue["C"].get_vector()
        
        psi = calc_dihedral(sN, sCA, sC, nN)
        omega = calc_dihedral(sCA, sC, nN, nCA)
        phi = calc_dihedral(sC, nN, nCA, nC)
        return phi, psi, omega
    
    def get_sine_cos_of_angels(self, pdb_file, residue_num, chain_id="A"):
        """Computes sin and cos of backbone dihedral angles: phi, psi and omega

        Args:
            pdb_file (str): file path
            residue_num (int): [description]
            chain_id (str, optional): [description]. Defaults to "A".

        Returns:
            float: sin and cos of phi, psi and omega
        """
        phi, psi, omega = self.get_angles(pdb_file=pdb_file, residue_num=residue_num, chain_id=chain_id)
        
        sin_phi, sin_psi, sin_omega = math.sin(phi), math.sin(psi), math.sin(omega)
        cos_phi, cos_psi, cos_omega = math.cos(phi), math.cos(psi), math.cos(omega)
        return sin_phi, sin_psi, sin_omega, cos_phi, cos_psi, cos_omega

# pdb_file = "data/pdbs_clean/1a5eA.pdb"        
# bd = BackboneDihedral()
# phi, psi, omega = bd.get_angles(pdb_file=pdb_file, residue_num=121, chain_id="A")
# print(phi, psi, omega)

# sin_phi, sin_psi, sin_omega, cos_phi, cos_psi, cos_omega = bd.get_sine_cos_of_angels(pdb_file=pdb_file, residue_num=121, chain_id="A")
# print(sin_phi, sin_psi, sin_omega, cos_phi, cos_psi, cos_omega)

