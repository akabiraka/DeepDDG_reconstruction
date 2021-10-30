import sys
sys.path.append("../DeepDDG_reconstruction")

import subprocess

class HydrogenBond(object):
    def __init__(self, output_dir=None) -> None:
        super().__init__()
        self.output_dir = "data/hbs/" if output_dir is None else output_dir
        self.hbplus_exe = "3rd_party_items/hbplus_exe"
        

    def set_up(self, pdb_file):  
        command = "cd {} && ../../{} ../../{}".format(self.output_dir, self.hbplus_exe, pdb_file)
        # print(command)
        subprocess.getoutput(command)

pdb_file = "data/pdbs_clean/1a43A.pdb" 
HB = HydrogenBond()
HB.set_up(pdb_file)