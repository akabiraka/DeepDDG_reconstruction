import sys
sys.path.append("../DeepDDG_reconstruction")

import pandas as pd
import numpy as np
from Bio.PDB import Polypeptide

from objects.PDBData import PDBData
from objects.Selector import ChainAndAminoAcidSelect
from biophysical_properties.BackboneDihedral import BackboneDihedral
from biophysical_properties.SASA import SASA
from biophysical_properties.SecondaryStructure import SecondaryStructure
from biophysical_properties.PSSM import PSSM

# configurations
pdb_dir = "data/pdbs/"
pdbs_clean_dir = "data/pdbs_clean/"
fastas_dir = "data/fastas/"
CIF = "mmCif"
input_file_path = "data/dataset_2.xlsx"
sheet_name = "train"
n_rows_to_skip = 0
n_rows_to_evalutate = 1

# object initialization
PDBData = PDBData(pdb_dir=pdb_dir)
backbone_dihedral = BackboneDihedral()
sasa = SASA()
secondary_structure = SecondaryStructure() 
pssm = PSSM()

# data generation
dfs = pd.read_excel(input_file_path, sheet_name=sheet_name)
# print(dfs.columns)

# helper utility functions, need to be moved later
def get_absolute_mutation_site(clean_pdb_file, mutation_site):
    """This is 0-based index. When a pdb's residue index starts from 
    not 1, use this. Need this when handling Fasta, SS type files."""
    starting_residue_index = int(open(clean_pdb_file).readline().split()[5])
    abs_mutation_site = mutation_site - starting_residue_index
    return abs_mutation_site

def create_mutant_fasta_file(wild_fasta_file, mutation_site, wild_residue):
    pdbid = wild_fasta_file.split("/")[2].split(".")[0]
    wild_residue = Polypeptide.three_to_one(wild_residue) if len(wild_residue)==3 else wild_residue
    with open(wild_fasta_file, "r") as wild_fasta_reader:
        lines = wild_fasta_reader.readlines()
        # print(lines)
        with open("{}{}_mutant.fasta".format(fastas_dir, pdbid), 'w') as mutant_fasta_writer:
            # print(lines[1][mutation_site])# = wild_residue
            fasta = lines[1][:mutation_site] + wild_residue + lines[1][mutation_site+1:]
            mutant_fasta_writer.write(lines[0].rstrip()+":mutant\n")
            mutant_fasta_writer.write(fasta)



for i, row in dfs.iterrows():
    if i+1 <= n_rows_to_skip: continue
    
    pdb_id = row["pdb_id"].lower()
    chain_id = "A"
    mutation_site = int(row["mutation_site"])
    wild_residue = row["wild_residue"]
    mutant_residue = row["mutant_residue"]

    PDBData.download_structure(pdb_id=pdb_id)
    clean_wild_protein_structure = PDBData.clean(pdb_id=pdb_id, chain_id=chain_id, selector=ChainAndAminoAcidSelect(chain_id))
    clean_pdb_file = pdbs_clean_dir+pdb_id+chain_id+".pdb"

    abs_mutation_site = get_absolute_mutation_site(clean_pdb_file, mutation_site)
    print("Protein no: {} -> {}, {}, {}".format(i+1, pdb_id, mutation_site, abs_mutation_site))

    PDBData.generate_fasta_from_pdb(pdb_id, chain_id, clean_pdb_file, save_as_fasta=True, output_fasta_dir="data/fastas/")
    wild_fasta_file = fastas_dir+pdb_id+chain_id+".fasta"
    create_mutant_fasta_file(wild_fasta_file, abs_mutation_site, mutant_residue)
    mutant_fasta_file = fastas_dir+pdb_id+chain_id+"_mutant.fasta"
    
    angles = backbone_dihedral.of_a_residue(pdb_file=clean_pdb_file, residue_num=mutation_site, chain_id=chain_id, type="both")

    sasa.set_up(pdb_file=clean_pdb_file)
    sasa_value = sasa.of_a_residue(residue_index=mutation_site)

    ss_one_hot = secondary_structure.of_a_residue(pdb_file=clean_pdb_file, residue_index=abs_mutation_site, type="one-hot")

    pssm.set_up(wild_fasta_file)
    wild_residue_pssm_score = pssm.of_a_residue(residue_index=abs_mutation_site)
    pssm.set_up(mutant_fasta_file)
    mutant_residue_pssm_score = pssm.of_a_residue(residue_index=abs_mutation_site)

    wild_residue_type = np.array(list('{0:05b}'.format(Polypeptide.three_to_index(wild_residue))), dtype=np.float32)
    mutant_residue_type = np.array(list('{0:05b}'.format(Polypeptide.three_to_index(mutant_residue))), dtype=np.float32)

    target_residue_features = np.concatenate((angles, sasa_value, ss_one_hot, wild_residue_pssm_score, mutant_residue_pssm_score, wild_residue_type, mutant_residue_type))

    print(target_residue_features.shape)
    print(angles, sasa_value, ss_one_hot, wild_residue_pssm_score, mutant_residue_pssm_score, wild_residue_type, mutant_residue_type)

    
    print()
    if i+1 == n_rows_to_skip+n_rows_to_evalutate: 
        break