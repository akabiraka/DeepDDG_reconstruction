import sys
sys.path.append("../DeepDDG_reconstruction")

from Bio.Blast.Applications import NcbipsiblastCommandline

class PSSM(object):
    """PSSSM stands for position-specific socring-matrix
    """
    
    def __init__(self) -> None:
        super().__init__()
        self.pssm_output_dir = "data/pssms/"
        self.psiblast_exe = "3rd_party_items/ncbi-blast-2.12.0+/bin/psiblast"
        self.db = "3rd_party_items/swissprot_db/swissprot"
        
        
    def compute(self, fasta_file):
        """
        This blast run the query sequence 3 iterations against a db using psiblast program,
        and save the output file in psiblast directory as xml format. 
        """
        # out_xml = CONFIGS.PSIBLAST_DIR + pdb_with_chain_id + CONFIGS.DOT_XML
        # out_pssm = CONFIGS.PSIBLAST_DIR + "pssm.txt"
        # out_ascii_pssm = CONFIGS.PSIBLAST_DIR + pdb_with_chain_id + ".pssm"
        pdbid = fasta_file.split("/")[2].split(".")[0]
        print("Computing PSSM for {} using psi-blast ... ...".format(pdbid))
        
        out_ascii_pssm = self.pssm_output_dir + pdbid +".pssm"
        E_VALUE_TRESH = "10"
        
        cline = NcbipsiblastCommandline(cmd=self.psiblast_exe, db=self.db, query=fasta_file,\
                                        evalue = E_VALUE_TRESH, outfmt = 5, num_iterations=3,\
                                        save_pssm_after_last_round=True, out_ascii_pssm=out_ascii_pssm)
                                        # out = out_xml, out_pssm=out_pssm, 
        cline()
        

fasta_file = "data/fastas/4eiuA.fasta"
pssm = PSSM()
pssm.compute(fasta_file)