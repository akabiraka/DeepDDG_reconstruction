import sys
sys.path.append("../DeepDDG_reconstruction")

import torch
from torch.utils.data import DataLoader
from DeepDDG.models import SRP, FCNNS
from DeepDDG.dataset import DeepDDGDataset

from analyzers.plotting import *

def test(data_loader):
    losses = []
    pred_ddgs = []
    exp_ddgs = []
    for i, data in enumerate(data_loader):
        pair_tensors, ddg = data
        exp_ddgs.append(ddg.item())
        pair_tensors = pair_tensors.to(device=device)
        ddg = ddg.to(device=device) / 10
        # print(pair_tensors.dtype, pair_tensors.shape)
        # print(ddg.dtype, ddg.shape)
        
        # running the model
        srp_outs = srp_model(pair_tensors)
        # print(srp_outs.shape) #batch_size,15,20
        ddg_pred = fcnn_model(srp_outs)
        # print(ddg_pred.shape) #batch_size,1
        pred_ddgs.append(ddg_pred.item()*10)
        # computing loss, backpropagate and optimizing model
        loss = criterion(ddg, ddg_pred)
        
        losses.append(loss.item())
        
    return losses, exp_ddgs, pred_ddgs


run_no = 1
N_neighbors=15

device = "cuda" if torch.cuda.is_available() else "cpu"
srp_model = SRP(in_features=51, out_features=20).to(device)
fcnn_model = FCNNS(in_features=N_neighbors*20).to(device)
criterion = torch.nn.MSELoss()

srp_model.load_state_dict(torch.load("outputs/models/run_1_srp_model_0.pth"))
fcnn_model.load_state_dict(torch.load("outputs/models/run_1_fcnn_model_0.pth"))
srp_model.eval()
fcnn_model.eval()

print("loading training dataset ... ...")
train_dataset = DeepDDGDataset(file="data/dataset_4_train_keep.csv",  data_dir="data/features_train/")
train_loader = DataLoader(train_dataset, batch_size=1, shuffle=False)
print("train dataset len:", train_dataset.__len__())
print("train loader size:", len(train_loader))
print("successfully loaded train dataset ... ...")

train_losses, train_exp_ddgs, train_pred_ddgs = test(train_loader)
plot_losses(train_losses, filename="run_{}_train_loss".format(run_no))
plot_pred_vs_exp_ddg(train_exp_ddgs, train_pred_ddgs, filename="run_1_train_pred_vs_exp_ddg".format(run_no))
plot_deviation_from_expected(train_exp_ddgs, train_pred_ddgs, filename="run_1_train_deviation_from_exp_ddg".format(run_no))
print("run_{}_train CC values".format(run_no))
compute_correlation_coefficient(train_exp_ddgs, train_pred_ddgs)


print("loading testing dataset ... ...")
test_dataset = DeepDDGDataset(file="data/dataset_4_test_keep.csv",  data_dir="data/features_test/")
test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)
print("test dataset len:", test_dataset.__len__())
print("test loader size:", len(test_loader))
print("successfully loaded testing dataset ... ...")

test_losses, test_exp_ddgs, test_pred_ddgs = test(test_loader)
plot_losses(test_losses, filename="run_{}_test_loss".format(run_no))
plot_pred_vs_exp_ddg(test_exp_ddgs, test_pred_ddgs, filename="run_1_test_pred_vs_exp_ddg".format(run_no))
plot_deviation_from_expected(test_exp_ddgs, test_pred_ddgs, filename="run_1_test_deviation_from_exp_ddg".format(run_no))
print("run_{}_test CC values".format(run_no))
compute_correlation_coefficient(test_exp_ddgs, test_pred_ddgs)