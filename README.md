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
