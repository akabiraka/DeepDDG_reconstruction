import sys
sys.path.append("../DeepDDG_reconstruction")

import torch
from torch import nn
import torch.optim as optim

class FCNNS(nn.Module):
    def __init__(self, in_features, dropout_probability=0.15) -> None:
        super().__init__()
        self.flatten = nn.Flatten()
        self.fcnns = nn.Sequential(
            nn.Linear(in_features, 100),
            nn.Dropout(p=dropout_probability),
            nn.ReLU(inplace=True),
            # nn.Linear(100, 100),
            # nn.Dropout(p=dropout_probability),
            # nn.ReLU(inplace=True),
            # nn.Linear(100, 100),
            # nn.Dropout(p=dropout_probability),
            # nn.ReLU(inplace=True),
            
        )
        self.regressor = nn.Sequential(
            nn.Linear(100, 1),
            nn.Dropout(p=dropout_probability)
        )
        self.classifier = nn.Sequential(
            nn.Linear(100, 1),
            nn.Dropout(p=dropout_probability)
        )
        self.softsign = nn.Softsign()

    def forward(self, x):
        # print(torch.sum(x, 1).shape)
        x = torch.sum(x, 1)
        # x = self.flatten(x)
        # print(x.shape)
        x = self.fcnns(x)
        # print(x.shape)
        x1 = self.regressor(x)
        pred_ddgs = self.softsign(x1)
        x2 = self.classifier(x)
        pred_classes = self.softsign(x2)
        return pred_ddgs, pred_classes

class SRP(nn.Module):
    """SRP for Shared Residue Pair network"""

    def __init__(self, in_features, out_features, dropout_probability=0.15) -> None:
        super().__init__()
        self.srp = nn.Sequential(
            nn.Linear(in_features, 100),
            nn.Dropout(p=dropout_probability),
            nn.ReLU(inplace=True),
            # nn.Linear(100, 100),
            # nn.Dropout(p=dropout_probability),
            # nn.ReLU(inplace=True),
            # nn.Linear(100, 100),
            # nn.Dropout(p=dropout_probability),
            # nn.ReLU(inplace=True),
            nn.Linear(100, out_features),
            nn.Dropout(p=dropout_probability),
            nn.ReLU(inplace=True),
        )

    def forward(self, x):
        # print(x.shape)
        x = self.srp(x)
        return x

def count_params(model):
    total_params_count = sum(p.numel() for p in model.parameters())
    trainable_params_count = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print("Total params: ", total_params_count, "Trainable params: ", trainable_params_count)

# device = 'cuda' if torch.cuda.is_available() else 'cpu'
# print('Using {} device'.format(device))

# # setting up the variables
# srp_model = SRP(in_features=2*25, out_features=20).to(device)
# count_params(srp_model)
# deepddg_model = DeepDDG(in_features=15*20).to(device)
# count_params(deepddg_model)
# criterion = nn.MSELoss()
# srp_optimizer = optim.Adam(srp_model.parameters(), lr=0.0008)
# deepddg_optimizer = optim.Adam(deepddg_model.parameters(), lr=0.0008)
# n_neighbors = 15

# # running the model
# srp_outs = []
# for n in range(n_neighbors):
#     target_neighbor_residue_features = torch.rand(1, 2, 25, device=device)
#     srp_outs.append(srp_model(target_neighbor_residue_features))

# concatenated_srp_outs = torch.cat(srp_outs, dim=1)
# # print(concatenated_srp_outs.shape)
# ddg_pred = deepddg_model(concatenated_srp_outs)
# # print(ddg_pred.shape, ddg_pred)

# # computing loss, backpropagate and optimizing model
# ddg_target = torch.randn(1, 1, device=device)
# loss = criterion(ddg_target, ddg_pred)
# print(loss)
# loss.backward()
# srp_optimizer.step()
# deepddg_optimizer.step()