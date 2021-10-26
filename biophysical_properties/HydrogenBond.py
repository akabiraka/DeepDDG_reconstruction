import sys
sys.path.append("../DeepDDG_reconstruction")

import subprocess

class HydrogenBond(object):
    def __init__(self) -> None:
        super().__init__()
        self.hbplus = "3rd_party_items/hbplus"
        

    def set_up(self, pdb_file, output_path):  
        command = "cd {} && ../../{} ../../{}".format(output_path, self.hbplus, pdb_file)
        # print(command)
        subprocess.getoutput(command)

pdb_file = "data/pdbs_clean/1a43A.pdb" 
output_path = "data/hbs/"
HB = HydrogenBond()
HB.set_up(pdb_file, output_path)