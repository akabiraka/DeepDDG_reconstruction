import sys
sys.path.append("../DeepDDG_reconstruction")


import requests
from Bio.PDB import *
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Polypeptide import PPBuilder 
from Bio import SeqIO

class PDBData(object):
    def __init__(self, pdb_dir="data/pdbs/"):
        super(PDBData, self).__init__()
        self.pdb_dir = pdb_dir
        self.pdb_format = "mmCif"
        self.pdb_ext = ".cif"
        self.parser = MMCIFParser(QUIET=True)
        self.pdbl = PDBList()
        self.pdbio = PDBIO()
        self.ppb=PPBuilder()
        
    def validate_pdb_id(self, pdb_id):
        """
        Validate the given pdb_id which must be a 4 letter pdb_id. 
        Else raise an Exception.
        """
        if len(pdb_id) == 4:
            return pdb_id
        elif len(pdb_id) > 4:
            pdb_id = pdb_id[0:4].lower()
        else:
            raise Exception("Reported ", pdb_id, "is not a pdb id in 4 letter format.")
        return pdb_id
        
    def download_structure(self, pdb_id):
        """Download the pdb structure into mmCif format.

        Args:
            pdb_id (string): pdb_id
        """
        self.pdbl.retrieve_pdb_file(pdb_id, pdir=self.pdb_dir, file_format=self.pdb_format)
    
    def get_first_chain_id(self, pdb_id):
        pdb_file = self.pdb_dir + pdb_id + self.pdb_ext
        structure = self.parser.get_structure("", pdb_file)[0]
        # print(list(structure.get_chains()))
        return list(structure.get_chains())[0].id
    
    def clean(self, pdb_id, chain_id, selector=None, clean_pdb_dir="data/pdbs_clean/"):
        """
            Given a Select instance, this function reads a mmCif protein data and saves into
            pdb format.
        Args:
            pdb_id ([string]): [description]
            chain_id ([character]): [description]
            selector: Subclass of Bio.PDB.Select class

        Returns:
            Structure: Return the cleaned structure
        """
        print("Cleaning {}:{} ... ..".format(pdb_id, chain_id))
        pdb_filename = self.pdb_dir + pdb_id + self.pdb_ext
        structure = self.parser.get_structure(pdb_id, pdb_filename)[0]
        self.pdbio.set_structure(structure)
        pdb_filename = clean_pdb_dir + pdb_id+chain_id + ".pdb"
        if selector is None:
            self.pdbio.save(pdb_filename)
        else:
            self.pdbio.save(pdb_filename, select=selector)
        return PDBParser(QUIET=True).get_structure(pdb_id, pdb_filename)
            
    def get_chain_from_structure(self, structure, chain_id):
        """Given a structure and chain_id, it returns 
        corresponding chain of the structure

        Args:
            structure ([Bio.PDB.Structure]): [An instance of Bio.PDB.Structure]
            chain_id ([character]): [chain_id]

        Returns:
            chain: Bio.PDB.Chain
        """
        models = list(structure.get_models())
        chains = list(models[0].get_chains())
        for chain in chains:
            if chain.id == chain_id:
                return chain
            
    def get_chain_from_clean_pdb(self, pdb_id, chain_id, pdb_filename):
        """Return chain given a pdb file.

        Args:
            pdb_id (string)
            chain_id (character)
            pdb_filename (path): A pdb file path. i.e: "data/pdbs_clean/1c8cA.pdb"
        Returns:
            Chain: A chain of a Protein data
        """
        structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_filename)
        models = list(structure.get_models())
        chains = list(models[0].get_chains())
        # for each chain
        for chain in chains:
            if chain.id == chain_id:
                return chain
            
    def generate_fasta_from_pdb(self, pdb_id, chain_id, input_pdb_filepath, save_as_fasta=False, output_fasta_dir=None):
        """Return sequence and length, and create fasta file from pdb data.

        Args:
            pdb_id (string)
            chain_id (string)
            input_pdb_filepath (string): the path must have a pdb file. i.e data/pdbs_clean/1amqA.pdb
        Returns:
            seq: primary sequence of the pdb file
        """
        
        structure = PDBParser(QUIET=True).get_structure(pdb_id, input_pdb_filepath)
        residues = structure[0][chain_id].get_residues()
        seq = ""
        for residue in residues:
            seq += Polypeptide.three_to_one(residue.get_resname())
        
        if save_as_fasta and output_fasta_dir is not None:
            with open("{}{}{}.fasta".format(output_fasta_dir, pdb_id, chain_id), "w") as fasta_file_handle:
                fasta_file_handle.write(">{}:{}\n".format(pdb_id.upper(), chain_id))
                fasta_file_handle.write(seq)
        print("Generating fasta {}:{}, Seq-len:{} ... ..".format(pdb_id, chain_id, len(seq)))
        return seq, len(seq)
    
    def create_mutant_fasta_file(self, wild_fasta_file, mutant_fasta_file, mutation_site, wild_residue):
        pdbid = wild_fasta_file.split("/")[2].split(".")[0]
        wild_residue = Polypeptide.three_to_one(wild_residue) if len(wild_residue)==3 else wild_residue
        with open(wild_fasta_file, "r") as wild_fasta_reader:
            lines = wild_fasta_reader.readlines()
            # print(lines)
            with open(mutant_fasta_file, 'w') as mutant_fasta_writer:
                # print(lines[1][mutation_site])# = wild_residue
                fasta = lines[1][:mutation_site] + wild_residue + lines[1][mutation_site+1:]
                mutant_fasta_writer.write(lines[0].rstrip()+":mutant\n")
                mutant_fasta_writer.write(fasta)
    
    def __save_fasta(self, pdb_id, fasta_text, is_save_file=True, fasta_dir="data/fastas/"):
        """Private method for download_fasta method

        Args:
            pdb_id (string): pdb_id
            fasta_text (string): text in fasta format
            is_save_file (bool, optional): [description]. Defaults to True.
            fasta_dir (str, optional): Directory path. Defaults to "data/fastas/".
        """
        if is_save_file:
            filepath = "{}{}_downloaded.fasta".format(fasta_dir, pdb_id)
            with open(filepath, "w") as fasta_file_handle:
                fasta_file_handle.write(fasta_text)
                fasta_file_handle.close()
        
    def download_fasta(self, pdb_id, is_save_file=True, fasta_dir="data/fastas/"):
        """Download fasta from pdb using pdb_id

        Args:
            pdb_id (string): pdb_id
            is_save_file (bool, optional): Self-explained. Defaults to True.
            fasta_dir (str, optional): Directory path. Defaults to "data/fastas/".
        """
        print("Downloading fasta {} ... ..".format(pdb_id))
        http = requests.Session()
        save_fasta_hook = lambda response, *args, **kwargs: self.__save_fasta(pdb_id, response.text, is_save_file, fasta_dir)
        http.hooks["response"] = save_fasta_hook
        http.get('https://www.rcsb.org/fasta/entry/'+pdb_id)
        
    def get_fasta_seq_record(self, pdb_id, fasta_dir="data/fastas/"):
        """Returns a Bio.SeqRecord.

        Args:
            pdb_id (string): pdb_id
            fasta_dir (string): directory path

        Returns:
            [Bio.SeqRecord]: returns the 1st Bio.SeqRecord from fasta file.
        """
        filepath = "{}{}_downloaded.fasta".format(fasta_dir, pdb_id)
        return next(SeqIO.parse(filepath, "fasta"))
        
    
    def get_seq_from_pdb(self, input_pdb_filepath):   
        """Given a input pdb filepath, it returns the sequence

        Args:
            input_pdb_filepath (string): example: "data/pdbs_clean/1amqA.pdb"

        Returns:
            Bio.Seq.Seq
        """
        structure = PDBParser(QUIET=True).get_structure("id", input_pdb_filepath)
        polypeptide = self.ppb.build_peptides(structure).__getitem__(0)
        return polypeptide.get_sequence()
    
    def get_first_residue_id(self, pdb_file, chain_id):
        residue_list = list(PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id].get_residues())
        _, residue_id, _ = residue_list[0].id
        return residue_id
                
    def get_last_residue_id(self, pdb_file, chain_id):
        residue_list = list(PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id].get_residues())
        _, residue_id, _ = residue_list[len(residue_list)-1].id
        return residue_id
    
    def get_residue_ids_dict(self, pdb_file, chain_id):
        residues = PDBParser(QUIET=True).get_structure("", pdb_file)[0][chain_id].get_residues()
        residue_ids_dict = {residue.id[1]:i for i, residue in enumerate(residues)}
        return residue_ids_dict
    
    
# clean_pdb_file = "data/pdbs_clean/1amqA.pdb" 
# PDBData = PDBData(pdb_dir="data/pdbs/")
# residue_ids_dict = PDBData.get_residue_ids_dict(pdb_file=clean_pdb_file, chain_id="A")
# print(residue_ids_dict)
# print(get_starting_residue_index(pdb_file=clean_pdb_file))

# print(get_last_residue_id("data/pdbs_clean/1a43A.pdb", "A"))
# print(get_first_residue_id("data/pdbs_clean/1a43A.pdb", "A"))