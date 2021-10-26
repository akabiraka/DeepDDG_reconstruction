# DeepDDG reconstruction

**[DeepDDG](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00697)**: Predicting the Stability Change of Protein Point Mutations Using Neural Networks

## Dataset

* *dataset_1:* is generated from *original_dataset_oct_2020* by changing the column names as "_" separated with lower case letters except for "pH" and "T". All the columns are kept as they were in the original one.
* *dataset_2*: the *mutation* column is seprated into three columns as *wild_residue*, *mutaion_site* and *mutant_residue*. *Wild* and *mutant* residue are 3-letter amino acid representation.
* The train set contains 5444 mutation data points from 211 proteins.
* The test set contains 276 mutation data points from 37 proteins.
* The train and test set contains no common proteins.
* The *output_images/data_distribution* directory compares the train and test set distributions upon:
  * Amino acid vs number of mutations
  * ddG vs number of mutations
  * Mutation site vs number of mutations
  * pH vs number of mutations
  * Proteins vs number of mutations
  * Temperature vs number of mutations

## Feature computation

* Backbone dihedral: phi, psi and omega angels of backbone atoms [N, CA, C] of a residue. (6 features)
  * Followed the below two links:
    * https://biopython.org/docs/dev/api/Bio.PDB.internal_coords.html#Bio.PDB.internal_coords.IC_Residue.pick_angle
    * https://biopython.org/docs/latest/api/Bio.PDB.vectors.html?highlight=calc_dihedral#Bio.PDB.vectors.calc_dihedral
* Solvent accessible surface area (SASA): using Naccess
  * http://www.bioinf.manchester.ac.uk/naccess/
  * https://www.biostars.org/p/43000/
* Secondary structures (SS): using Stride
  * Downloaded and compiled following http://webclu.bio.wzw.tum.de/stride/
  * DeepDDG mentions they use 3 bits for SS, but did not mention for which 3 types.
  * So, I keep coil (C), beta-sheet (E) and alpha-helix (H). Although DSSP classifies SS as 8 types.
* Hydrogen bonds: using HBPLUS
  * To download: https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/
  * Usage manual: https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html
  *
* Position-specific scoring matrix (PSSM): using PSI-BLAST (authors version 2.7.1)
  * Download the latest version ncbi-blast executables from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
    * File name: ncbi-blast-2.12.0+-x64-linux.tar.gz
    * ncbi-blast executables contain blastp, blastn, psiblast and so on.
    * Quick start on BLAST: https://www.ncbi.nlm.nih.gov/books/NBK569856/
  * Download the database from https://ftp.ncbi.nlm.nih.gov/blast/db/
    * File name: swissprot.tar.gz (because it is small to setup and test)
  * Authors' used database is *rp-seq-55*
    * Will be setup later
    *
  * Authors' also applied softmax over PSSM row-wise as the final PSSM.
  * Biopython provides `NcbipsiblastCommandline` to use psi-blast software.
    * `from Bio.Blast.Applications import NcbipsiblastCommandline`
  * Sample example of BLAST usage from command line:
    * To download swissprot: `3rd_party_items/ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress swissprot`
    * To download all "nr" databases: "nr" is the non-redundent version of the sequence database.
      * `3rd_party_items/ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress nr [*]`
    * To create blast-db from fasta sequences:
      * `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in path/4eiuA.fasta -input_type fasta -out path_to_save/db_name`
* Pairwise fitness score (PFS): using HH-suite3.0
  * HH-suite homepage: https://github.com/soedinglab/hh-suite
  * To download precompiled HH-suite: https://mmseqs.com/hhsuite/
    * File name: hhsuite-linux-sse2.tar.gz
  * User guide: https://github.com/soedinglab/hh-suite/wiki#generating-a-multiple-sequence-alignment-using-hhblits
  * To set up the database I first download the scop40_01Mar17.tgz from https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/ (because it is small 381M)
    * Sample usage: ./3rd_party_items/hhsuite/bin/hhblits -i data/fastas/4eiuA.fasta -o test.msa -d ./3rd_party_items/scop40_01Mar17/scop40
  * To download uniprot_20: https://wwwuser.gwdg.de/~compbiol/data/hhsuite/benchmark/
* Multiple sequence alignment (MSA):
*

1. 1A43 -> if I want to take a chain, I shall take the "A"
2. 1SPB:P -> I thought P is a chain but not!!! what is this then?
3. 1A7V:Q1A -> the Q at 1st residue will be replace by A, need confirmation.
4. 1ACB:I:F10W -> there is no I chain,
5. 1CEY_F14N -> this is not same as 1A7V:Q1A since 14th residue is not F
6. 1CEY_WT
7. 1CFD_1-75
8. 1CFD_1-78_F19Y
9. 1VII_N68AK70M
