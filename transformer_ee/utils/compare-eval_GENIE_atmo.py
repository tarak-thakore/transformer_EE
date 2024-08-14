import os
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy import stats
from scipy.optimize import curve_fit
import binstat



# Check if variable exists (optional)
# if 'model_path' in os.environ:  # Replace with the variable name you want to check
#     model_path = os.environ['model_path']
#     print("model_path : ",model_path)
# else:
#     print("model_path not found")

# filepath = model_path + "/result.npz"
model1_name = "GENIEv3-0-6-Honda-Truth-hA-LFG_wLeptonScalars_EnergyOnly_MAE"
model2_name = "GENIEv3-0-6-Honda-Truth-hA-LFG_wLeptonScalars_MAE"
file1path = "/home/jbarrow/MLProject2/save/model/GENIEv3-0-6-Honda-Truth-hA-LFG_wLeptonScalars_EnergyOnly_MAE/model_GENIEv3-0-6-Honda-Truth-hA-LFG_wLeptonScalars_EnergyOnly_MAE/result.npz"
file2path = "/home/jbarrow/MLProject2/save/model/GENIEv3-0-6-Honda-Truth-hA-LFG_wLeptonScalars_MAE/model_GENIEv3-0-6-Honda-Truth-hA-LFG_wLeptonScalars_MAE/result.npz"
plotpath = "/home/jbarrow/MLProject2/save/model/Comparisons/"

print("Contents of the npz file1:")
with np.load(file1path) as file1:
  for key in file1.keys():
      print(key)

print("Contents of the npz file2:")
with np.load(file2path) as file2:
  for key in file2.keys():
      print(key) 
      
file1 = np.load(file1path)
trueval_EOnly = file1['trueval']
prediction_EOnly = file1['prediction']

file2 = np.load(file2path)
trueval_All = file2['trueval']
prediction_All = file2['prediction']

print("trueval_EOnly shape: ",trueval_EOnly.shape,prediction_EOnly.shape)
n_val1,dim1 = trueval_EOnly.shape

print("trueval_All shape: ",trueval_All.shape,prediction_All.shape)
n_val2,dim2 = trueval_All.shape

#True Variables
E_nu_true_EOnly = trueval_EOnly[:,0]
E_nu_true_All = trueval_All[:,0]
# Px_nu_true_EOnly = trueval_EOnly[:,1]
# Py_nu_true_EOnly = trueval_EOnly[:,2]
# Pz_nu_true_EOnly = trueval_EOnly[:,3]
# #True angle variables
# Cos_Theta_nu_true_EOnly = Py_nu_true_EOnly / (((Px_nu_true_EOnly) ** 2 + (Py_nu_true_EOnly) ** 2 + (Pz_nu_true_EOnly) **2) ** 0.5 )
# #Theta_nu_true_EOnly = math.acos(Cos_Theta_nu_true_EOnly)

#Predicted variables
E_nu_pred_EOnly = prediction_EOnly[:,0]
E_nu_pred_All = prediction_All[:,0]
# Px_nu_pred_EOnly = prediction_EOnly[:,1]
# Py_nu_pred_EOnly = prediction_EOnly[:,2]
# Pz_nu_pred_EOnly = prediction_EOnly[:,3]
# #Predicted angle variables
# Cos_Theta_nu_pred_EOnly = Py_nu_pred_EOnly / (((Px_nu_pred_EOnly) ** 2 + (Py_nu_pred_EOnly) ** 2 + (Pz_nu_pred_EOnly) **2) ** 0.5 )
# #Theta_nu_pred_EOnly = math.acos(Cos_Theta_nu_pred_EOnly)

binstat.plot_xstat(E_nu_true_EOnly,E_nu_pred_EOnly,name=plotpath +"/en",title=r'Transformer Performance for $E_{\nu}$ Reconstruction',scale='linear',xlabel=r'True $E_{\nu}$',ylabel=r'Predicted $E_{\nu}$')
binstat.plot_xstat(E_nu_true_EOnly,(E_nu_pred_EOnly-E_nu_true_EOnly)/E_nu_true_EOnly,name=plotpath +"/en",title="Atmospheric Neutrino Energy Reconstruction",scale='linear',xlabel="True Neutrino Energy",ylabel="Resolution on Neutrino Energy")

#binstat.plot_xstat(Cos_Theta_nu_true_EOnly,Cos_Theta_nu_pred_EOnly,name=plotpath +"/ct",title="Transformer Performance for $cos(\theta_{\nu})$",xlabel="True $cos(\theta_{\nu})$",ylabel="Predicted $cos(\theta_{\nu})$")
#binstat.plot_xstat(Cos_Theta_nu_true_EOnly,Cos_Theta_nu_pred_EOnly,name=plotpath +"/ct",title="Atmospheric Neutrinos' Cosine of Incoming Angle",xlabel="True Cosine of Incoming Angle",ylabel="Predicted Cosine of Incoming Angle")

binstat.plot_y_hist(100.*((E_nu_pred_EOnly-E_nu_true_EOnly)/E_nu_true_EOnly),100.*((E_nu_pred_All-E_nu_true_All)/E_nu_true_All),name=plotpath+"/en_res",labels=["Energy Only", "Composite"],xrange=(-100.0,100.0),yrange=(0.5,1.0E6),xlabel="Percent Difference Energy Resolution",log=True,vline=True,bins=200,colors=["red","blue"])

#binstat.plot_2d_hist_count(E_nu_true_EOnly,E_nu_pred_EOnly,name=plotpath +"/en_hist2d",xrange=(0.0,4.0),yrange=(0.0,4.0),title="Transformer Performance for $E_{\nu}$ Reconstruction",scale='linear',xlabel="True $E_{\nu}$",ylabel="Predicted $E_{\nu}$")
#binstat.plot_2d_hist_count(E_nu_true_EOnly,E_nu_pred_EOnly,name=plotpath +"/en_hist2d",xrange=(0.0,4.0),yrange=(0.0,4.0),title="Atmospheric Neutrino Energy Reconstruction",scale='linear',zscale='log',xlabel="True Neutrino Energy",ylabel="Predicted Neutrino Energy")

#binstat.plot_2d_hist_count(Px_nu_true_EOnly,Px_nu_pred_EOnly,name=model_name +"/ct_hist2d",xrange=(-1.,1.),yrange=(-1.,1.))

#binstat.plot_2d_hist_count(Cos_Theta_nu_true_EOnly,Cos_Theta_nu_pred_EOnly,name=plotpath +"/ct_hist2d",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Transformer Performance for $cos(\theta_{\nu})$ Reconstruction",xlabel="True cos(\theta_{\nu})",ylabel="Predicted cos(\theta_{\nu})")
#binstat.plot_2d_hist_count(Cos_Theta_nu_true_EOnly,Cos_Theta_nu_pred_EOnly,name=plotpath +"/ct_hist2d",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Atmospheric Neutrinos' Cosine of Incoming Angle",scale='linear',zscale='log',xlabel="True Cosine of Incoming Angle",ylabel="Predicted Cosine of Incoming Angle")