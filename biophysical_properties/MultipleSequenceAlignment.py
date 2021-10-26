import sys
sys.path.append("../DeepDDG_reconstruction")
import os
import subprocess

class MultipleSequenceAlignment(object):
    def __init__(self, db=None) -> None:
        super().__init__()
        self.msa_output_dir = "data/msas/"
        self.hhblits_exe = "./3rd_party_items/hhsuite/bin/hhblits"
        self.db = "./3rd_party_items/scop40_01Mar17/scop40" if db is None else db
        
        
    def set_db(self, db):
        self.db = db
        
    def set_up(self, fasta_file, force=False):
        pdbid = fasta_file.split("/")[2].split(".")[0]
        output_file_path = self.msa_output_dir + pdbid +".msa"
        if os.path.exists(output_file_path) and force==False: 
            print("MSA is already set up for {}. To set-up again, set force=True".format(pdbid))
            return
        else:
            print("Computing MSA for {} using HH-suite ... ...".format(pdbid))
            command = "{} -i {} -o {} -d {}".format(self.hhblits_exe, fasta_file, output_file_path, self.db)
            subprocess.getoutput(command)
        
# 
fasta_file = "data/fastas/4eiuA.fasta"
MSA = MultipleSequenceAlignment()
# MSA.set_up(fasta_file)
# MSA.set_up(fasta_file, force=True)