import sys
sys.path.append("../DeepDDG_reconstruction")

import subprocess
import os

class SASA(object):
    """SASA stands for Solvent Accessible Surface Area.
    """
    def __init__(self, output_dir=None) -> None:
        super().__init__()
        self.output_dir = "data/sasas/" if output_dir is None else output_dir
        self.naccess_exe = "3rd_party_items/naccess"

    def set_up(self, pdb_file, force=False):  
        pdbid = pdb_file.split("/")[2].split(".")[0]

        if os.path.exists(self.output_dir+pdbid+".asa") and force==False: 
            print("SASA is already set up for {}. To set-up again, set force=True".format(pdbid))
            return
        else:
            print("Computing SASA for {} using psi-blast ... ...".format(pdbid)) 
            command = "cd {} && ../../{} ../../{}".format(self.output_dir, self.naccess_exe, pdb_file)
            subprocess.getoutput(command)
            subprocess.getoutput("cd {} && mv .asa {}.asa".format(self.output_dir, pdbid))
            subprocess.getoutput("cd {} && rm .rsa && rm .log".format(self.output_dir))


pdb_file = "data/pdbs_clean/1a43A.pdb" 
sasa = SASA()
sasa.set_up(pdb_file)
# sasa.set_up(pdb_file, force=True)