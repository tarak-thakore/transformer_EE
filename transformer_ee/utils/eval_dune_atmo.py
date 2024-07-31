import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import binstat



# Check if variable exists (optional)
if 'model_path' in os.environ:  # Replace with the variable name you want to check
    model_path = os.environ['model_path']
    print("model_path : ",model_path)
else:
    print("model_path not found")

filepath = model_path + "/result.npz"

print("Contents of the npz file:")
with np.load(filepath) as file:
  for key in file.keys():
      print(key)  
      
file = np.load(filepath)
trueval = file['trueval']
prediction = file['prediction']

print("trueval shape: ",trueval.shape,prediction.shape)
n_val,dim = trueval.shape

en_true = trueval[:,0]
ct_true = trueval[:,1]
en_pred = prediction[:,0]
ct_pred = prediction[:,1]

resolution_en = (en_pred - en_true)/en_true
resolution_ct = (ct_pred - en_true)/ct_true

binstat.plot_xstat(en_true,resolution_en,name=model_path +"/en",title=r"$E_\nu$",scale='linear',xlabel=r"$E_{true}$",ylabel=r"$E_{resolution}$")
binstat.plot_xstat(ct_true,resolution_ct,name=model_path +"/ct",title=r"$\theta_\nu$ (degrees)",xlabel=r"$\theta_{true}$ (degrees)",ylabel=r"$\theta_{resolution}$ (degrees)")
binstat.plot_y_hist(resolution_en,name=model_path+"/en_res")
#binstat.plot_y_hist(resolution_ct,name=model_path+"/ct_res")
binstat.plot_2d_hist_count(en_true,en_pred,name=model_path +"/en_hist2d",xrange=(0,10),yrange=(0,10),title=r"$E_\nu$",scale='linear',xlabel=r"$E_{true}$",ylabel=r"$E_{rec}$")
#binstat.plot_2d_hist_count(ct_true,ct_pred,name=model_path +"/ct_hist2d",xrange=(-1.,1.),yrange=(-1.,1.))
binstat.plot_2d_hist_count(ct_true,ct_pred,name=model_path +"/ct_hist2d",xrange=(0.,180.),yrange=(0.,180.),title=r"$\theta_\nu$ (degrees) ",xlabel=r"$\theta_{true} (degrees)$",ylabel=r"$\theta_{rec}$ (degrees)")
