import sys
sys.path.append("../DeepDDG_reconstruction")

import subprocess
import os
import pandas as pd
import numpy as np

class SASA(object):
    """SASA stands for Solvent Accessible Surface Area.
    """
    def __init__(self, output_dir=None) -> None:
        super().__init__()
        self.output_dir = "data/sasas/" if output_dir is None else output_dir
        self.naccess_exe = "3rd_party_items/naccess"
        self.pdb_id = None

    def set_up(self, pdb_file, force=False):  
        pdbid = pdb_file.split("/")[2].split(".")[0]
        self.pdb_id = pdbid

        if os.path.exists(self.output_dir+pdbid+".rsa") and force==False: 
            print("SASA is already set up for {}. To set-up again, set force=True.".format(pdbid))
            return
        else:
            print("Computing SASA for {} using Naccess ... ...".format(pdbid)) 
            command = "cd {} && ../../{} ../../{}".format(self.output_dir, self.naccess_exe, pdb_file)
            subprocess.getoutput(command)
            subprocess.getoutput("cd {} && mv .rsa {}.rsa".format(self.output_dir, pdbid))
            subprocess.getoutput("cd {} && rm .asa && rm .log".format(self.output_dir))

    def parse_rsa_file(self, rsa_file):
        df = pd.read_csv(rsa_file, delim_whitespace=True, header=None, skiprows=[0, 1, 2, 3]) # skipping 1st 4 rows
        df = df.head(-4) # removing last 4 rows
        df[3] = df[3].astype(int)
        # df[3] = df[3]-df[3][0]+1
        # print(df)
        return df

    def of_a_residue(self, residue_index):
        rsa_file = self.output_dir + self.pdb_id + ".rsa"
        asa_df = self.parse_rsa_file(rsa_file)
        asa_value = asa_df[asa_df[3]==residue_index][4].astype(float)
        # print(np.array(asa_value))
        return np.array(asa_value, dtype=np.float32)

# pdb_file = "data/pdbs_clean/1a43A.pdb" 
# sasa = SASA()
# sasa.set_up(pdb_file)
# sasa_value = sasa.of_a_residue(residue_index=1)
# print(sasa_value)