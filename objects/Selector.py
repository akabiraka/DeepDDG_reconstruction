__author__ = "Anowarul Kabir"
__updated__ = "2020-08-29 20:13:28"

import sys
sys.path.append("../DeepDDG_reconstruction")

from Bio.PDB import *
from Bio.PDB.PDBIO import Select

class CAAtomSelector(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    
    def  accept_chain(self, chain):
        # print(chain.id)
        if chain.id == self.chain_id:
            return 1
        else:
            return 0
        
    def accept_atom(self, atom):
        """Overload this to reject atoms for output."""
        if atom.name == "CA":
            return 1
        else:
            return 0        
        
class StandardAminoAcidSelector(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
        
    def  accept_chain(self, chain):
        # print(chain.id, self.chain_id)
        if self.chain_id == chain.id:
            return 1
        else:
            return 0
            
    def accept_residue(self, residue):
        if residue.get_resname() in standard_aa_names:
            return 1
        else:
            return 0

class AllBackboneAtomSelector(Select):
    """Backbone atoms: CA, CB, N, O

    Args:
        Select ([type]): [description]
    """
    def __init__(self, chain_id):
        self.chain_id = chain_id
    
    def  accept_chain(self, chain):
        # print(chain.id)
        if chain.id == self.chain_id:
            return 1
        else:
            return 0
        
    def accept_atom(self, atom):
        """Overload this to reject atoms for output."""
        if atom.name == "CA":
            return 1
        elif atom.name == "CB":
            return 1
        elif atom.name == "N":
            return 1
        elif atom.name == "O":
            return 1
        else:
            return 0
        
class ChainAndAminoAcidSelect(Select):
    def __init__(self, chain_id):
        super(ChainAndAminoAcidSelect).__init__()
        self.chain_id = chain_id
        
    def  accept_chain(self, chain):
        if chain.id == self.chain_id:
            # print(chain.id, self.chain_id)
            return 1
        else:
            return 0
            
    def accept_residue(self, residue):
        _, residue_id, insertion_code = residue.id
        if residue.get_resname() in standard_aa_names and insertion_code==" ":
            # print(residue.get_resname())
            return 1
        else:
            return 0