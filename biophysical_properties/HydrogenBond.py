import sys
sys.path.append("../DeepDDG_reconstruction")

import subprocess
import numpy as np
import os

class HydrogenBond(object):
    def __init__(self, output_dir=None) -> None:
        super().__init__()
        self.output_dir = "data/hbs/" if output_dir is None else output_dir
        self.hbplus_exe = "3rd_party_items/hbplus"
        self.pdb_id = None
        

    def set_up(self, pdb_file, force=False):  
        self.pdb_id = pdb_file.split("/")[2].split(".")[0]

        if os.path.exists(self.output_dir+self.pdb_id+".hb2") and force==False: 
            print("Hydrogen bond is already set up for {}. To set-up again, set force=True.".format(self.pdb_id))
            return
        else:
            command = "cd {} && ../../{} ../../{}".format(self.output_dir, self.hbplus_exe, pdb_file)
            # print(command)
            subprocess.getoutput(command)

    def get(self, target_residue_id, neighbor_residue_id):
        hb_file = self.output_dir+self.pdb_id+".hb2"
        num_of_hydrogen_bonds = np.array([0, 0, 0, 0], dtype=np.float32)
        with open(hb_file, "r") as hb_file_reader:
            lines = hb_file_reader.readlines()[8:]
            for line in lines:
                donar_residue_id = int(line[2:5])
                acceptor_residue_id = int(line[15:19])
                donar_atom_category = line[33]
                acceptor_atom_category = line[34]

                if donar_residue_id==neighbor_residue_id and acceptor_residue_id==target_residue_id:
                    if donar_atom_category=="M" and acceptor_atom_category=="M": num_of_hydrogen_bonds[0] += 1
                    if donar_atom_category=="M" and acceptor_atom_category=="S": num_of_hydrogen_bonds[1] += 1
                    if donar_atom_category=="S" and acceptor_atom_category=="M": num_of_hydrogen_bonds[2] += 1
                    if donar_atom_category=="S" and acceptor_atom_category=="S": num_of_hydrogen_bonds[3] += 1


        return num_of_hydrogen_bonds


# pdb_file = "data/pdbs_clean/1a43A.pdb" 
# target_residue_id=184
# neighbor_residue_id=183

# HB = HydrogenBond()
# HB.set_up(pdb_file)
# num_of_hydrogen_bonds = HB.get(target_residue_id, neighbor_residue_id)
# print(num_of_hydrogen_bonds)