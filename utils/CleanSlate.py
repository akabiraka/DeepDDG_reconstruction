__author__ = "Anowarul Kabir"
__updated__ = "2020-07-24 22:16:53"

import sys
sys.path.append("../protein_data")

import os

class CleanSlate(object):
    def __init__(self):
        super(CleanSlate, self).__init__()
        
    def clean_all(self):
        data_root = "data/"
        self.clean_all_files(mydir=data_root + "fastas/", ext=".fasta")
        self.clean_all_files(mydir=data_root + "features/", ext=".pt")
        self.clean_all_files(mydir=data_root + "hbs/", ext=".hb2")
        self.clean_all_files(mydir=data_root + "msas/", ext=".a3m")
        self.clean_all_files(mydir=data_root + "pdbs/", ext=".cif")
        self.clean_all_files(mydir=data_root + "pdbs_clean/", ext=".pdb")
        self.clean_all_files(mydir=data_root + "pssms/", ext=".pssm")
        self.clean_all_files(mydir=data_root + "sasas/", ext=".rsa")
        self.clean_all_files(mydir=data_root + "secondary_structures/", ext=".ss")
        
    def delete_file(self, filepath):
        if os.path.exists(filepath):
            os.remove(filepath)

    def clean_all_files(self, mydir, ext):
        """
        dir: "/data/pdbs"
        ext: ".cif"
        """
        for f in os.listdir(mydir): 
            if f.endswith(ext): 
                os.remove(os.path.join(mydir, f))


# cln = CleanSlate()
# cln.clean_all()