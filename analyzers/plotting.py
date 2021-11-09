import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr

def plot_losses(train_losses, validation_losses, filename):    
    plt.plot(train_losses, c="lightgreen", label="Train loss")
    plt.plot(validation_losses, c="salmon", label="Validation loss")
    plt.xlabel("Number of epochs")
    plt.ylabel("Mean-squared error (MSE)")
    plt.savefig("outputs/images/model_analysis/{}.pdf".format(filename), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.cla()
    # plt.show()

def plot_pred_vs_exp_ddg(exp_ddgs, pred_ddgs, filename):
    x = np.array(exp_ddgs)
    y = np.array(pred_ddgs)
    plt.scatter(x, y, c="salmon", marker=".")

    m, b = np.polyfit(x, y, 1) #m = slope, b=intercept
    plt.plot(x, m*x + b, c="lightgreen")

    plt.xlabel("Expected ddG (kcal/mol)")
    plt.ylabel("Predicted ddG (kcal/mol)")
    plt.savefig("outputs/images/model_analysis/{}.pdf".format(filename), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.cla()
    # plt.show()
    
def plot_deviation_from_expected(exp_ddgs, pred_ddgs, filename):
    x = np.array(exp_ddgs)
    y = np.array(pred_ddgs)
    error = x-y
    
    zeros = np.zeros_like(x)
    plt.errorbar(x, zeros, yerr=error, fmt='o', color='green',
             ecolor='salmon', elinewidth=3, capsize=0, alpha=.5)

    plt.xlabel("Expected ddG (kcal/mol)")
    plt.ylabel("Error deviation (kcal/mol)")
    plt.savefig("outputs/images/model_analysis/{}.pdf".format(filename), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.cla()
    # plt.show()

def compute_correlation_coefficient(exp_ddgs, pred_ddgs):    
    pcc_r, pcc_p_value = pearsonr(exp_ddgs, pred_ddgs)
    scc_r, scc_p_value = spearmanr(exp_ddgs, pred_ddgs)       
    print("PCC: ", pcc_r, pcc_p_value)
    print("SCC: ", scc_r, scc_p_value) 
    
# plot_losses(train_losses, filename="run_1_train_loss")
# plot_pred_vs_exp_ddg(train_exp_ddgs, train_pred_ddgs, filename="run_1_train_pred_vs_exp_ddg")
# plot_deviation_from_expected(train_exp_ddgs, train_pred_ddgs, filename="run_1_train_deviation_from_exp_ddg")
# print("run_1_train CC values")
# compute_correlation_coefficient(train_exp_ddgs, train_pred_ddgs)