import pandas as pd

train_set_df = pd.read_excel("data/dataset_1.xlsx", sheet_name="train")
test_set_df = pd.read_excel("data/dataset_1.xlsx", sheet_name="test")

train_set_df = train_set_df[["pdb_id", "mutation", "ddG"]] 
grouped_train_set_df = train_set_df.groupby(by=["pdb_id", "mutation"], as_index=False).mean()
# print(grouped_train_set_df)
grouped_train_set_df.to_excel("data/dataset_3.xlsx", sheet_name="train", index=False)


test_set_df = test_set_df[["pdb_id", "mutation", "ddG"]]
grouped_test_set_df = test_set_df.groupby(by=["pdb_id", "mutation"], as_index=False).mean()
# print(grouped_test_set_df)
with pd.ExcelWriter('data/dataset_3.xlsx', mode="a") as writer: 
    grouped_test_set_df.to_excel(writer, sheet_name="test", index=False)

unique_proteins_train = grouped_train_set_df["pdb_id"].drop_duplicates().to_list()
unique_proteins_test = grouped_test_set_df["pdb_id"].drop_duplicates().to_list()
print(len(unique_proteins_train), len(unique_proteins_test)) # 209, 37

common_proteins = list(set(unique_proteins_train) & set(unique_proteins_test))
print(len(common_proteins)) # 0