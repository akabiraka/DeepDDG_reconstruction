import sys
sys.path.append("../DeepDDG_reconstruction")

import numpy as np
import torch
from torch import nn, optim
from torch.utils.data import DataLoader
from DeepDDG.models import SRP, FCNNS, count_params
from DeepDDG.dataset import DeepDDGDataset
import matplotlib.pyplot as plt
from analyzers.plotting import *
from datasets.validation_set_generator import *

def validate():
    srp_model.eval()
    fcnn_model.eval()
    losses = []
    for i, data in enumerate(validation_loader):
        pair_tensors, ddG, exp_classes = data
        pair_tensors = pair_tensors.to(device=device)
        ddG = ddG.to(device=device)
        srp_outs = srp_model(pair_tensors)
        ddg_pred, pred_classes = fcnn_model(srp_outs)
        mse_loss = mse_criterion(ddG, ddg_pred*10)
        l1_loss = l1_criterion(exp_classes, pred_classes)
        loss = alpha*mse_loss+beta*l1_loss
        losses.append(loss)
    return torch.stack(losses).mean().item()

def train():
    srp_model.train()
    fcnn_model.train()
    losses = []
    for i, data in enumerate(train_loader):
        pair_tensors, ddG, exp_classes = data
        pair_tensors = pair_tensors.to(device=device)
        ddG = ddG.to(device=device)
        # print(pair_tensors.dtype, pair_tensors.shape)
        # print(ddG.dtype, ddG.shape)
        # print(class_label.dtype, class_label.shape)
        
        
        # running the model
        srp_model.zero_grad()
        fcnn_model.zero_grad()
        srp_outs = srp_model(pair_tensors)
        ddg_pred, pred_classes = fcnn_model(srp_outs)
        # print(srp_outs.shape) #batch_size,15,20
        # print(ddg_pred.shape) #batch_size,1
        
        # computing loss, backpropagate and optimizing model
        mse_loss = mse_criterion(ddG, ddg_pred*10)
        l1_loss = l1_criterion(exp_classes, pred_classes)
        loss = alpha*mse_loss+beta*l1_loss
        
        # print(loss)
        loss.backward()
        optimizer.step()
        losses.append(loss)
        
    return torch.stack(losses).mean().item()

run_no = 11
learning_rates = [0.001]
for i, learning_rate in enumerate(learning_rates):
    print("ith_start=", i)
    print("initializing variables ... ...")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    batch_size = 200
    n_epochs = 50
    N_neighbors = 15
    # eps = 0.001
    weight_decay = 0.0008
    alpha = 1
    beta = 0.01
    print("run_no=", run_no)
    print("batch_size=", batch_size)
    print("n_epochs=", n_epochs)
    print("lr=", learning_rate) 
    print("N_neighbors=", N_neighbors)
    # print("eps=", eps)
    print("weight_decay=", weight_decay)
    
    print("initializing models, loss function and optimizers ... ...") 
    srp_model = SRP(in_features=51, out_features=20, dropout_probability=0.5).to(device)
    fcnn_model = FCNNS(in_features=20, dropout_probability=0.5).to(device)
    print(srp_model)
    print(fcnn_model)
    mse_criterion = nn.MSELoss()
    l1_criterion = nn.L1Loss()
    optimizer = optim.Adam([{"params":srp_model.parameters()}, {"params": fcnn_model.parameters()}], lr=learning_rate, weight_decay=weight_decay)
    
    print("Printing the number of parameters of SRP and FCNN.")
    count_params(srp_model)
    count_params(fcnn_model)
    
    # generate_validation_set(n_pdbids=20)
    
    print("loading training dataset ... ...")
    train_dataset = DeepDDGDataset(file="data/dataset_5_train.csv", data_dir="data/features_train/")
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    print("train dataset len:", train_dataset.__len__())
    print("train loader size:", len(train_loader))
    print("successfully loaded training dataset ... ...")
    
    print("loading validation dataset ... ...")
    validation_dataset = DeepDDGDataset(file="data/dataset_5_validation.csv", data_dir="data/features_train/")
    validation_loader = DataLoader(validation_dataset, batch_size=validation_dataset.__len__(), shuffle=True)
    print("validation dataset len:", validation_dataset.__len__())
    print("validation loader size:", len(validation_loader))
    print("successfully loaded validation dataset ... ...")
    
    print("training models ... ...")
    best_loss = np.inf        
    train_losses = []
    validation_losses = []
    for epoch in range(1, n_epochs+1):
        train_loss = train()    
        train_losses.append(train_loss)
        validation_loss = validate()
        validation_losses.append(validation_loss)
        print("[{}/{}] train-loss: {:.4f}, validation-loss: {:.4f}".format(epoch, n_epochs, train_loss, validation_loss))
            
        # if epoch%10==0:
        #     plt.plot(train_losses)
        #     plt.plot(validation_losses)
        #     plt.show()
            
        if validation_loss < best_loss:
            torch.save(srp_model.state_dict(), "outputs/models/run_{}_srp_model_{}.pth".format(run_no, i))
            torch.save(fcnn_model.state_dict(), "outputs/models/run_{}_fcnn_model_{}.pth".format(run_no, i))
        
        # print(optimizer.param_groups[0]['lr'])            
    # print("train_losses_by_epochs=", train_losses)
    plot_losses(train_losses, validation_losses, filename="run_{}_train_loss".format(run_no))