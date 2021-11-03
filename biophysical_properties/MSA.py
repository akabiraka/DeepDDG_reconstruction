import sys
sys.path.append("../DeepDDG_reconstruction")
import os
import subprocess

class MSA(object):
    """MSA for Multiple Sequence Alignment
    """
    
    def __init__(self, db=None, output_dir=None) -> None:
        super().__init__()
        self.db = "3rd_party_items/scop40_01Mar17/scop40" if db is None else db
        self.output_dir = "data/msas/" if output_dir is None else output_dir
        self.hhblits_exe = "3rd_party_items/hhsuite/bin/hhblits"
        
        
    def set_up(self, fasta_file, force=False):
        pdbid = fasta_file.split("/")[2].split(".")[0]
        
        # for raw outputs
        # output_file_path = self.output_dir + pdbid +".msa"
        # command = "{} -i {} -o {} -d {}".format(self.hhblits_exe, fasta_file, output_file_path, self.db)
        
        # for MSA only
        output_file_path = self.output_dir + pdbid +".a3m"
        command = "{} -i {} -oa3m {} -d {}".format(self.hhblits_exe, fasta_file, output_file_path, self.db)
        
        if os.path.exists(output_file_path) and force==False: 
            print("MSA is already set up for {}. To set-up again, set force=True".format(pdbid))
            return
        else:
            print("Computing MSA for {} using HH-suite ... ...".format(pdbid))
            subprocess.getoutput(command)
            # subprocess.getoutput("rm data/fastas/{}.hhr".format(pdbid))
        
# -oa3m query.a3m
# fasta_file = "data/fastas/4eiuA.fasta"
# msa = MSA(db="3rd_party_items/uniprot20_2015_06_bench/uniprot20_2015_06/uniprot20_2015_06_merged/uniprot20_2015_06_merged")
# msa.set_up(fasta_file)
# # MSA.set_up(fasta_file, force=False)