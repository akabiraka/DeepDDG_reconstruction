import sys
sys.path.append("../DeepDDG_reconstruction")

import torch
from torch.utils.data import DataLoader
from DeepDDG.models import SRP, FCNNS
from DeepDDG.dataset import DeepDDGDataset

import numpy as np
from analyzers.plotting import *

def test(data_loader):
    mse_losses = []
    mae_losses = []
    pred_ddgs = []
    exp_ddgs = []
    for i, data in enumerate(data_loader):
        pair_tensors, ddg = data
        exp_ddgs.append(ddg.item())
        pair_tensors = pair_tensors.to(device=device)
        ddg = ddg.to(device=device)
        # print(pair_tensors.dtype, pair_tensors.shape)
        # print(ddg.dtype, ddg.shape)
        
        # running the model
        srp_outs = srp_model(pair_tensors)
        # print(srp_outs.shape) #batch_size,15,20
        ddg_pred = fcnn_model(srp_outs)*10
        # print(ddg_pred.shape) #batch_size,1
        pred_ddgs.append(ddg_pred.item())
        # computing loss, backpropagate and optimizing model
        
        mse_loss = mse_criterion(ddg, ddg_pred)
        mae_loss = mae_criterion(ddg, ddg_pred)
        
        mse_losses.append(mse_loss)
        mae_losses.append(mae_loss)
    
    mse = torch.stack(mse_losses).mean().item()
    mae = torch.stack(mae_losses).mean().item()
    return mse, mae, exp_ddgs, pred_ddgs


run_no = 8
N_neighbors=15

device = "cuda" if torch.cuda.is_available() else "cpu"
srp_model = SRP(in_features=51, out_features=20).to(device)
fcnn_model = FCNNS(in_features=N_neighbors*20).to(device)
mse_criterion = torch.nn.MSELoss()
mae_criterion = torch.nn.L1Loss()

srp_model.load_state_dict(torch.load("outputs/models/run_{}_srp_model_0.pth".format(run_no)))
fcnn_model.load_state_dict(torch.load("outputs/models/run_{}_fcnn_model_0.pth".format(run_no)))
srp_model.eval()
fcnn_model.eval()

print("loading training dataset ... ...")
train_dataset = DeepDDGDataset(file="data/dataset_5_train.csv", data_dir="data/features_train/")
train_loader = DataLoader(train_dataset, batch_size=1, shuffle=True)
print("train dataset len:", train_dataset.__len__())
print("train loader size:", len(train_loader))
print("successfully loaded training dataset ... ...")

print("loading validation dataset ... ...")
validation_dataset = DeepDDGDataset(file="data/dataset_5_validation.csv", data_dir="data/features_train/")
validation_loader = DataLoader(validation_dataset, batch_size=1, shuffle=True)
print("validation dataset len:", validation_dataset.__len__())
print("validation loader size:", len(validation_loader))
print("successfully loaded validation dataset ... ...")

print("loading testing dataset ... ...")
test_dataset = DeepDDGDataset(file="data/dataset_4_test_keep.csv",  data_dir="data/features_test/")
test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)
print("test dataset len:", test_dataset.__len__())
print("test loader size:", len(test_loader))
print("successfully loaded testing dataset ... ...")

train_mse, train_mae, train_exp_ddgs, train_pred_ddgs = test(train_loader)
plot_pred_vs_exp_ddg(train_exp_ddgs, train_pred_ddgs, filename="run_{}_train_pred_vs_exp_ddg".format(run_no))
plot_deviation_from_expected(train_exp_ddgs, train_pred_ddgs, filename="run_{}_train_deviation_from_exp_ddg".format(run_no))

print("train MSE: ", np.mean(train_mse))
print("train MAE: ", np.mean(train_mae))
print("run_{}_train CC values".format(run_no))
compute_correlation_coefficient(train_exp_ddgs, train_pred_ddgs)

val_mse, val_mae, val_exp_ddgs, val_pred_ddgs = test(validation_loader)
plot_pred_vs_exp_ddg(val_exp_ddgs, val_pred_ddgs, filename="run_{}_val_pred_vs_exp_ddg".format(run_no))
plot_deviation_from_expected(val_exp_ddgs, val_pred_ddgs, filename="run_{}_val_deviation_from_exp_ddg".format(run_no))

print("val MSE: ", np.mean(val_mse))
print("val MAE: ", np.mean(val_mae))
print("run_{}_val CC values".format(run_no))
compute_correlation_coefficient(val_exp_ddgs, val_pred_ddgs)


test_mse, test_mae, test_exp_ddgs, test_pred_ddgs = test(test_loader)
plot_pred_vs_exp_ddg(test_exp_ddgs, test_pred_ddgs, filename="run_{}_test_pred_vs_exp_ddg".format(run_no))
plot_deviation_from_expected(test_exp_ddgs, test_pred_ddgs, filename="run_{}_test_deviation_from_exp_ddg".format(run_no))

print("test MSE: ", np.mean(test_mse))
print("test MAE: ", np.mean(test_mae))
print("run_{}_test CC values".format(run_no))
compute_correlation_coefficient(test_exp_ddgs, test_pred_ddgs)