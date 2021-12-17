# DeepDDG reconstruction

**[DeepDDG](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00697)**: Predicting the Stability Change of Protein Point Mutations Using Neural Networks

Further analysis of model building and evaluation result can be found at the [report](https://github.com/akabiraka/DeepDDG_reconstruction/blob/main/report/DeepDDG.pdf).

## Data cleaning

* *dataset_1:* is generated from *original_dataset_oct_2020* by changing the column names as "_" separated with lower case letters except for "pH" and "T". All the columns are kept as they were in the original one.
* *dataset_2*: the *mutation* column is seprated into three columns as *wild_residue*, *mutaion_site* and *mutant_residue*. *Wild* and *mutant* residue are 3-letter amino acid representation.
* The train set contains 5444 mutation data points from 209 proteins.
* The test set contains 276 mutation data points from 37 proteins.
* Different ddG values are reported for same mutation, which is removed in *dataset_3*.
* Now the train and test set contains 4344 mutations of 209 proteins and 253 mutations of 37 proteins
* The train and test set contains no common proteins.
* The *output_images/data_distribution* directory compares the train and test set distributions upon:
  * Amino acid vs number of mutations
  * ddG vs number of mutations
  * Mutation site vs number of mutations
  * pH vs number of mutations
  * Proteins vs number of mutations
  * Temperature vs number of mutations
* Some data issues are reported at the bottom of this file.

## Feature computation

* Backbone dihedral: phi, psi and omega angels of backbone atoms [N, CA, C] of a residue. (6 features)

  * Followed the below two links:
    * https://biopython.org/docs/dev/api/Bio.PDB.internal_coords.html#Bio.PDB.internal_coords.IC_Residue.pick_angle
    * https://biopython.org/docs/latest/api/Bio.PDB.vectors.html?highlight=calc_dihedral#Bio.PDB.vectors.calc_dihedral
* Solvent accessible surface area (SASA): using Naccess

  * To download: http://www.bioinf.manchester.ac.uk/naccess/
  * Get the access key from the authors, and compile it from /home directory. If the directory path contains space (i.e. ..../New volme/....) it cannot it compiled. It requires gfortran and csh.
  * Sample usage: naccess 1a43A.pdb
  * After installation the path should be: `"3rd_party_items/naccess"`
* Secondary structures (SS): using Stride

  * Downloaded and compiled following http://webclu.bio.wzw.tum.de/stride/
  * DeepDDG mentions they use 3 bits for SS, but did not mention for which 3 types.
  * So, I keep coil (C), beta-sheet (E) and alpha-helix (H). Although DSSP classifies SS as 8 types.
  * After installation the path should be: `"3rd_party_items/stride"`
* Hydrogen bonds: using HBPLUS

  * To download: https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/
  * Usage manual: https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/manual.html
  * After installation the path should be: `"3rd_party_items/hbplus"`
* Position-specific scoring matrix (PSSM): using PSI-BLAST (authors version 2.7.1)

  * Download the latest version ncbi-blast executables from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

    * File name: ncbi-blast-2.12.0+-x64-linux.tar.gz
    * ncbi-blast executables contain blastp, blastn, psiblast and so on.
    * Quick start on BLAST: https://www.ncbi.nlm.nih.gov/books/NBK569856/
    * After installation the path should be: `"3rd_party_items/ncbi-blast-2.12.0+/bin/psiblast"`
  * Database setup

    * I used rp-seq-15 instead of rp-seq-55, because rp-seq-55 is more than 20GB of size where as rp-seq-15 is less than 6GB  of size.
    * Download link: https://proteininformationresource.org/rps/
    * Make the database compatible for psi-blast: `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in path/to/rp-seq-55.fasta -input_type fasta -out 3rd_party_items/rp_seq_15/`
  * Since PSSM feature generation takes longer, this is generated distributedly in GMU argo cluster:

    * `sbatch jobs/distributed_pssm_generator.sh`
  * Authors' also applied softmax over PSSM row-wise as the final PSSM.
  * The following items are optional for practicing purposes:

    * Download the database from https://ftp.ncbi.nlm.nih.gov/blast/db/*

      * File name: swissprot.tar.gz (because it is small to setup and test)
    * Biopython provides `NcbipsiblastCommandline` to use psi-blast software.

      * `from Bio.Blast.Applications import NcbipsiblastCommandline`
    * Sample example of BLAST usage from command line:

      * To download swissprot: `3rd_party_items/ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress swissprot`
      * To download all "nr" databases: "nr" is the non-redundent version of the sequence database.
        * `3rd_party_items/ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress nr [*]`
        * `ncbi-blast-2.12.0+/bin/update_blastdb.pl --decompress nr.02`
      * To create blast-db from fasta sequences:
        * `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in path/4eiuA.fasta -input_type fasta -out path_to_save/db_name`
        * `3rd_party_items/ncbi-blast-2.12.0+/bin/makeblastdb -dbtype prot -in 3rd_party_items/rp-seqs-15.fasta -input_type fasta -out 3rd_party_items/rp_req_15/rp_req_15`
* Pairwise fitness score (PFS) after Multiple sequence alignment (MSA): using HH-suite3.0

  * HH-suite homepage: https://github.com/soedinglab/hh-suite
  * To download precompiled HH-suite: https://mmseqs.com/hhsuite/
    * File name: hhsuite-linux-sse2.tar.gz
  * User guide: https://github.com/soedinglab/hh-suite/wiki#generating-a-multiple-sequence-alignment-using-hhblits
  * To set up the database I first download the scop40_01Mar17.tgz from https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/ (because it is small 381M)
    * Sample usage: ./3rd_party_items/hhsuite/bin/hhblits -i data/fastas/4eiuA.fasta -o test.msa -d ./3rd_party_items/scop40_01Mar17/scop40
  * To download uniprot_20: https://wwwuser.gwdg.de/~compbiol/data/hhsuite/benchmark/
  * After installation the path should be: `"3rd_party_items/hhsuite/bin/hhblits"`
* After the avobe setup successfully, run the following command

  * `python datasets/input_data_generator.py`
  * To generate train/test set specific feature, line 19 and 20 must be changed accordingly.

## Training and testing

* To train the model, run: `python DeepDDG/train.py`
* To test the model, run: `python DeepDDG/test.py`

## Clarifications

* Column name: PDB ID with modifications to be made
* 1A43 -> if I want to take a chain, I shall take the "A"
* 1SPB:P -> I thought P is a chain but not!!! what is this then?
* 1A7V:Q1A -> the 1st residue is Q, which will be substituted by A.
* 1ACB:I:F10W -> there is no I chain,
* 1CEY_F14N -> you can see that although it's the 13th residue according to the PDB, the author labeled it 14.
* 1CEY_WT
* 1CFD_1-75: I think I need to consider only 1-75th residue
* 1CFD_1-78_F19Y: I think I need to consider 1-78th residue, here 19th is F, will be substituted by Y.
* 1VII_N68AK70M: N 68 ->A and K 70-> M
* The SASA value is taken from rsa file the absolute (ABS) all atoms value.

## Data issues

* For the same mutation of the same protein the dataset contains different ddg values. i.e 1A43 has 9 mutations where 4 of them are repeated with different ddG values. Took the average following Potapov et al (https://doi.org/10.1093/protein/gzp030).
* No chain id is not given in the original dataset. I selected the 1st chain by default. Because in many cases chain A is not present. And in many cases 1st chain is not A. Such cases are:
  * ```
    PDBID Chain-id
    1lmb     4
    1tup     A
    1azp     A
    1bf4     A
    1otr     B
    1glu     A
    1hcq     A
    1iv7     B
    ```
* 1hfzA has the 1st residue as 1x. I take the starting index as where the residue id becomes integer.
* Entry 2A01 was removed. Check this: https://www.rcsb.org/structure/removed/2A01
* Last residue does not have any dihedral angles. Therefore, last residue cannot be neighbor residue.
* 1a7cA does not have residue from 334 to 347.
* Some ddG values are out of range ([-10, 10]).
* 1am7A does not have residue of number 17.
* 1amqA does not have residue of number 407.
* In many cases, a particular neighbor residue does not have Ca, C or N atoms which is required to compute dihedral angles. To solve this problem, this neighbor residue is avoided and set the next residue as neighbor and follow the process again.
* 2ptlA_T_19_A: this does not have any hit while computing PSSM, need bigger database, and no PSSM file generated. Therefore returned 0 for softmax pssm.
* 4hxj_A_D_141_A: the secondary structure is not of the same size of the number of residues, manually corrected.
* 5np8_A_T_378_P: this does not have residue from 373 to 381. But 369th residue is Thr (T). So manually changed mutation point.
* 2arf_A_H_1069_Q: when SASA is computed there is no gap between chain_id and residue_num column. This is corrected manually.
* Note that: this type of errors may occur in other cases, which is not reported since the adapted code successfully avoided or solved them.
