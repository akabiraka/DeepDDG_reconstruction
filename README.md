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

* Backbone dihedral: sin and cos of phi, psi and omega angels of backbone atoms [N, CA, C] of a residue. (6 features)
* Solvent accessible surface area (SASA):
* Secondary structures (SS):
* Hydrogen bonds:
* Position-specific scoring matrix (PSSM):
* Pairwise fitness score (PFS):
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
