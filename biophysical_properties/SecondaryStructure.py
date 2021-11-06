import sys
sys.path.append("../DeepDDG_reconstruction")

import os
import subprocess
import string
import numpy as np

class SecondaryStructure(object):
    def __init__(self, output_dir=None) -> None:
        super().__init__()
        self.output_dir = "data/secondary_structures/" if output_dir is None else output_dir
        self.stride_exe = "3rd_party_items/stride"
        self.SS_dictionary = "CEH"
        
    def __parse_stride_file(self, file):
        f = open(file, "r")
        start_flag = False
        ss = ""
        for line in f:
            if start_flag:
                ss=ss+line.split()[5]
            if line.__contains__("|--Structure--|"): start_flag = True
        return ss
                
    def __get_by_Stride(self, pdb_file, save=False):
        pdb_file_name = pdb_file.split("/")[2].split(".")[0]
        output_file_path = self.output_dir + pdb_file_name + ".ss"
        result = subprocess.run([self.stride_exe, pdb_file, "-f"+output_file_path], stdout=subprocess.PIPE)
        ss = self.__parse_stride_file(output_file_path)
        if save:
            with open(output_file_path, "w") as f:
                f.write(">"+pdb_file_name+"\n")
                f.write(ss)
        else: os.remove(output_file_path)
        return ss
    
    def save(self, pdb_file):
        return self.__get_by_Stride(self, pdb_file, save=True)
    
    def set_secondary_structure_dictionary(self, ss_dict):
        self.SS_dictionary = ss_dict
    
    def of_a_residue(self, pdb_file, residue_index, return_type="one-hot"):
        """Returns secondary structure of a residue.

        Args:
            pdb_file (str): a pdb file path
            residue_index (int): It's 0 based. 
            return_type (str): Could be "one-hot" or "letter"
        Returns:
            str: a letter defined in self.SS_dictionary
        """
        ss = self.__get_by_Stride(pdb_file, save=True)
        letter = ss[residue_index]
        if return_type=="letter": return letter
        elif return_type=="one-hot": return np.array([0 if char != letter else 1 for char in self.SS_dictionary], dtype=np.float32)
    
    def of_some_residues(self, pdb_file, from_residue=0, n_residues=None, return_type="one-hot"):
        """Returns secondary structure of a set of residues defined by from_residue and n_residues.
        When n_residues is None, it returns the whole secondary strucutre of the protein as letters 
        defined by self.SS_dictionary.

        Args:
            pdb_file (str): a pdb file path.
            from_residue (int, optional): index of a residue. Defaults to 0.
            n_residues (int, optional): number of residues to process starting from from_residue. 
                Defaults to None.
            return_type (str): Could be "one-hot" or "letter"
        Returns:
            str: letters defined in self.SS_dictionary
        """
        ss = self.__get_by_Stride(pdb_file, save=False)
        ss = ss if n_residues is None else ss[from_residue: from_residue+n_residues-1]
        if return_type=="letter": return ss
        elif return_type=="one-hot":
            return np.array([[0 if char != letter else 1 for char in self.SS_dictionary] for letter in ss])

        
# clean_pdb_file = "data/pdbs_clean/1amqA.pdb" 
# secondaryStructure = SecondaryStructure() 
# result = secondaryStructure.of_a_residue(pdb_file=clean_pdb_file, residue_index=401-5-1, return_type="one-hot") # from opposite access check
# print(result)
# result = secondaryStructure.of_a_residue(pdb_file, 0, return_type="one-hot") # boundary value check
# result = secondaryStructure.of_a_residue(pdb_file, 155, return_type="one-hot") # boundary value check
# result = secondaryStructure.of_a_residue(pdb_file, 156, return_type="one-hot") # index out of range check
# result = secondaryStructure.of_a_residue(pdb_file, 5000, return_type="one-hot") # index out of range check

# result = secondaryStructure.of_a_residue(pdb_file, -1, return_type="letter") # from opposite access check
# result = secondaryStructure.of_a_residue(pdb_file, 0, return_type="letter") # boundary value check
# result = secondaryStructure.of_a_residue(pdb_file, 155, return_type="letter") # boundary value check
# result = secondaryStructure.of_a_residue(pdb_file, 156, return_type="letter") # index out of range check
# result = secondaryStructure.of_a_residue(pdb_file, 5000, return_type="letter") # index out of range check

# result = secondaryStructure.of_some_residues(pdb_file, return_type="one-hot") # getting all SS values check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=-10, n_residues=1, return_type="one-hot") # from opposite access check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=0, n_residues=10, return_type="one-hot") # regular check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=150, n_residues=6, return_type="one-hot") # boundary value check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=150, n_residues=500, return_type="one-hot") # index out of range check

# result = secondaryStructure.of_some_residues(pdb_file, return_type="letter") # getting all SS values check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=-10, n_residues=1, return_type="letter") # from opposite access check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=0, n_residues=10, return_type="letter") # regular check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=150, n_residues=6, return_type="letter") # boundary value check
# result = secondaryStructure.of_some_residues(pdb_file, from_residue=150, n_residues=500, return_type="letter") # index out of range check

# print(result)