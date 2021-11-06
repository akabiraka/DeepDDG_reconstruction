import sys
sys.path.append("../DeepDDG_reconstruction")

import math
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral
from utils import pdb_utils

class BackboneDihedral(object):
    def __init__(self) -> None:
        super().__init__()
        
    def of_a_residue(self, pdb_file, residue_num, chain_id="A", return_type=None):
        """
        Returns phi, psi and omega angles of a residue.
        # used the following two links to compute phi, psi and omega angels
            # https://biopython.org/docs/dev/api/Bio.PDB.internal_coords.html#Bio.PDB.internal_coords.IC_Residue.pick_angle
            # https://biopython.org/docs/latest/api/Bio.PDB.vectors.html?highlight=calc_dihedral#Bio.PDB.vectors.calc_dihedral
            # (sN, sCA, sC, nN)   # psi
            # (sCA, sC, nN, nCA)  # omega i+1
            # (sC, nN, nCA, nC)   # phi i+1

        Args:
            pdb_file ([return_type]): [description]
            residue_num ([return_type]): [description]
            chain_id (str, optional): [description]. Defaults to "A".
            return_type (str, optional): How the angles will be returned. 
                return_type can be "both", "sin", "cos" or None. 
                The default is phi, psi and omega values in a array.
        Returns:
            1D np array: angles.
        """
        last_residue_id = pdb_utils.get_last_residue_id(pdb_file, chain_id)
        # print(last_residue_id)
        if residue_num==last_residue_id:
            # since last residue does not have any dyhedral angles
            return np.array([0, 0, 0, 0, 0, 0], dtype=np.float32)
        
        chain = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id]
        residue = chain[residue_num]
        
        # finding next residue, since sometimes it is not possible to find the next residue by just next residue id.
        i=1
        while(1):
            if chain.has_id(residue_num+i):
                next_residue = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id][residue_num+i]
                break
            i+=1
        
        
        sN, sCA, sC = residue["N"].get_vector(), residue["CA"].get_vector(), residue["C"].get_vector()
        nN, nCA, nC = next_residue["N"].get_vector(), next_residue["CA"].get_vector(), next_residue["C"].get_vector()
        
        psi = calc_dihedral(sN, sCA, sC, nN)
        omega = calc_dihedral(sCA, sC, nN, nCA)
        phi = calc_dihedral(sC, nN, nCA, nC)
        
        if return_type=="sin": return np.array([math.sin(phi), math.sin(psi), math.sin(omega)], dtype=np.float32)
        elif return_type=="cos": return np.array([math.cos(phi), math.cos(psi), math.cos(omega)], dtype=np.float32)
        elif return_type=="both": return np.array([math.sin(phi), math.sin(psi), math.sin(omega), math.cos(phi), math.cos(psi), math.cos(omega)], dtype=np.float32)
        else: return np.array([phi, psi, omega], dtype=np.float32)
    
    
    def of_some_residues(self, pdb_file, chain_id="A", from_residue=1, n_residues=None, return_type=None):
        """It conputes angles for a subset of residues using from_residue and n_residues parameter.
        If they are not passed, the angles of all residues of the protein will be returned.

        Args:
            pdb_file (str): a pdb file path
            chain_id (str, optional): chain id. Defaults to "A".
            from_residue (int, optional): a residue number. Defaults to 1.
            n_residues (int, optional): how many residues to consider. Defaults to None.
            return_type (str, optional): how the angles will be returned. 
                How the angles will be returned. 
                return_type can be "both", "sin", "cos" or None. 

        Returns:
            2D np array: angles.
        """
        chain = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id]
        if n_residues is None: n_residues = len(list(chain.get_residues()))
        else: n_residues = from_residue+n_residues 
        all_angles = []
        for i in range(from_residue, n_residues):
            all_angles.append(self.of_a_residue(pdb_file, i, chain_id, return_type))
        
        return np.array(all_angles)
            
        
# pdb_file = "data/pdbs_clean/1am7A.pdb"        
# bd = BackboneDihedral()

# angles = bd.of_a_residue(pdb_file=pdb_file, residue_num=16, chain_id="A", return_type="both")
# print(angles)

# angles = bd.of_some_residues(pdb_file=pdb_file, chain_id="A", from_residue=150, n_residues=6, return_type="both")
# print(angles)
