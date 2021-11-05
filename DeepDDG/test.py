import sys
sys.path.append("../DeepDDG_reconstruction")

import torch
from torch.utils.data import DataLoader
from DeepDDG.models import SRP, FCNNS
from DeepDDG.dataset import DeepDDGDataset

N_neighbors=15

device = "cuda" if torch.cuda.is_available() else "cpu"
srp_model = SRP(in_features=51, out_features=20).to(device)
fcnn_model = FCNNS(in_features=N_neighbors*20).to(device)
criterion = torch.nn.MSELoss()

srp_model.load_state_dict(torch.load("outputs/models/srp_model_0.pth"))
fcnn_model.load_state_dict(torch.load("outputs/models/fcnn_model_0.pth"))
srp_model.eval()
fcnn_model.eval()


print("loading testing dataset ... ...")
test_dataset = DeepDDGDataset(file="data/dataset_4_train.csv", device=device)
test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)
print("test dataset len:", test_dataset.__len__())
print("test loader size:", len(test_loader))
print("successfully loaded testing dataset ... ...")


losses = []
for i, data in enumerate(test_loader):
    pair_tensors, ddG = data
    pair_tensors = pair_tensors.to(device=device)
    ddG = ddG.to(device=device) / 10
    # print(pair_tensors.dtype, pair_tensors.shape)
    # print(ddG.dtype, ddG.shape)
    
    # running the model
    srp_outs = srp_model(pair_tensors)
    # print(srp_outs.shape) #batch_size,15,20
    ddg_pred = fcnn_model(srp_outs)
    # print(ddg_pred.shape) #batch_size,1
    
    # computing loss, backpropagate and optimizing model
    loss = criterion(ddG, ddg_pred)
    
    losses.append(loss)
    print(loss)

print(torch.stack(losses).mean().item())
print("test_losses=", losses)