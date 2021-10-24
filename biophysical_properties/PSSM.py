import sys
sys.path.append("../DeepDDG_reconstruction")

import pandas as pd
from Bio.Blast.Applications import NcbipsiblastCommandline

class PSSM(object):
    """PSSSM stands for position-specific socring-matrix
    """
    def __init__(self, db=None) -> None:
        super().__init__()
        self.pssm_output_dir = "data/pssms/"
        self.psiblast_exe = "3rd_party_items/ncbi-blast-2.12.0+/bin/psiblast"
        self.db = "3rd_party_items/swissprot_db/swissprot" if db is None else db
        
        
    def set_db(self, db):
        self.db = db
        
        
    def set_up(self, fasta_file):
        """This blast run the query sequence 3 iterations against a db using psiblast program,
        and save the output file in pssms directory. 

        Args:
            fasta_file (str): file path
        """
        pdbid = fasta_file.split("/")[2].split(".")[0]
        print("Computing PSSM for {} using psi-blast ... ...".format(pdbid))
        
        output_file_path = self.pssm_output_dir + pdbid +".pssm"
        E_VALUE_TRESH = "10"
        
        cline = NcbipsiblastCommandline(cmd=self.psiblast_exe, db=self.db, query=fasta_file,\
                                        evalue = E_VALUE_TRESH, outfmt = 5, num_iterations=3,\
                                        save_pssm_after_last_round=True, out_ascii_pssm=output_file_path)# 
                                        # out = out_xml, out_pssm=out_pssm, out_ascii_pssm=output_file_path)
        cline()
        
    
    def __parse_pssm_output_file(self, pssm_file):
        df = pd.read_csv(pssm_file, delim_whitespace=True, header=None, skiprows=[1, 2]) # skipping 1st 2 rows
        df = df.head(-5) # removing last 5 rows
        return df
    
    
    def __typed_output(self, result, type):
        if type=="df": return result
        elif type=="np": return result.to_numpy()
    
    
    def of_a_residue(self, pssm_file, residue_index, type="np"):
        df = self.__parse_pssm_output_file(pssm_file)
        pssm = df.loc[residue_index, 2:21]
        return self.__typed_output(pssm, type)
    
    
    def of_some_residues(self, pssm_file, from_residue=0, n_residues=None, type="np"):
        df = self.__parse_pssm_output_file(pssm_file)
        pssm = df.loc[:, 2:21] if n_residues is None else df.loc[from_residue:from_residue+n_residues-1, 2:21]
        return self.__typed_output(pssm, type)
        

# fasta_file = "data/fastas/4eiuA.fasta"
# pssm_file = "data/pssms/4eiuA.pssm"

# sample usage
# pssm = PSSM()
# pssm.set_up(fasta_file)

# result = pssm.of_a_residue(pssm_file, 0, type="np") # check boundary value
# result = pssm.of_a_residue(pssm_file, 241, type="np") # check boundary value

# result = pssm.of_some_residues(pssm_file, from_residue=0, n_residues=None, type="np") # getting pssm for all residues
# result = pssm.of_some_residues(pssm_file, from_residue=0, n_residues=5, type="np")
# result = pssm.of_some_residues(pssm_file, from_residue=235, n_residues=100, type="np") # getting pssm upto last residue
# print(result.shape)