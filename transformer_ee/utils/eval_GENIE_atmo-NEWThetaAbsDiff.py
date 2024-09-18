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
#model_name = "GENIEv3-0-6-Honda-Truth-hA-LFG_wLeptonScalars_MAE"
filepath = "/exp/dune/app/users/rrichi/save/model/LossVars_Nu_Energy_Nu_Mom_X_Nu_Mom_Y_Nu_Mom_Z/GENIEv3-0-6-Honda-Truth-hA-LFG_AnyNu_Inclusive_Thresh_p1to1_eventnum_All_NpNpi_MSE/model_GENIEv3-0-6-Honda-Truth-hA-LFG_AnyNu_Inclusive_Thresh_p1to1_eventnum_All_NpNpi_MSE/result.npz"
plotpath = "/home/jbarrow/MLProject2/save/model/NEW_FullyVectorizedFinalStates/RichiWork_Vectorized/new_plots/LossVars_Nu_Energy_Nu_Mom_X_Nu_Mom_Y_Nu_Mom_Z/GENIEv3-0-6-Honda-Truth-hA-LFG_AnyNu_Inclusive_Thresh_p1to1_eventnum_All_NpNpi_MSE"

print("Contents of the npz file:")
with np.load(filepath) as file:
  for key in file.keys():
      print(key)  
      
file = np.load(filepath)
trueval = file['trueval']
prediction = file['prediction']

# print("trueval shape: ",trueval.shape,prediction.shape)
# n_val,dim = trueval.shape

print("Original trueval shape: ", trueval.shape)
print("Original prediction shape: ", prediction.shape)

# Apply the condition E_nu_true < 1.0
mask = trueval[:, 0] < 1.0
trueval = trueval[mask]
prediction = prediction[mask]

print("Filtered trueval shape: ", trueval.shape)
print("Filtered prediction shape: ", prediction.shape)
n_val, dim = trueval.shape

#True Variables
E_nu_true = trueval[:,0]
Px_nu_true = trueval[:,1]
Py_nu_true = trueval[:,2]
Pz_nu_true = trueval[:,3]
#True angle variables
Cos_Theta_nu_true = Py_nu_true / (((Px_nu_true) ** 2 + (Py_nu_true) ** 2 + (Pz_nu_true) **2) ** 0.5 )
Theta_nu_true = (180.0/np.pi)*np.arccos(Cos_Theta_nu_true)
#Theta_nu_true = math.arccos(Cos_Theta_nu_true)

#Predicted variables
E_nu_pred = prediction[:,0]
Px_nu_pred = prediction[:,1]
Py_nu_pred = prediction[:,2]
Pz_nu_pred = prediction[:,3]
#Predicted angle variables
Cos_Theta_nu_pred = Py_nu_pred / (((Px_nu_pred) ** 2 + (Py_nu_pred) ** 2 + (Pz_nu_pred) **2) ** 0.5 )
Theta_nu_pred = (180.0/np.pi)*np.arccos(Cos_Theta_nu_pred)
#Theta_nu_pred = math.arccos(Cos_Theta_nu_pred)

#Make resolution variables
#Energy
E_nu_res_percent=100.*((E_nu_pred-E_nu_true)/E_nu_true)
E_nu_res_percent_diff=100.*((E_nu_pred-E_nu_true)/((E_nu_true+E_nu_pred)/2))
E_nu_res_abs_diff=E_nu_pred-E_nu_true
#Theta
Theta_nu_res_percent=100.*((Theta_nu_pred-Theta_nu_true)/Theta_nu_true)
Theta_nu_res_percent_diff=100.*((Theta_nu_pred-Theta_nu_true)/((Theta_nu_true+Theta_nu_pred)/2))
Theta_nu_res_abs_diff=Theta_nu_pred-Theta_nu_true
#CosTheta
Cos_Theta_nu_res_percent=100.*((Cos_Theta_nu_pred-Cos_Theta_nu_true)/Cos_Theta_nu_true)
Cos_Theta_nu_res_percent_diff=100.*((Cos_Theta_nu_pred-Cos_Theta_nu_true)/((Cos_Theta_nu_true+Cos_Theta_nu_pred)/2))
Cos_Theta_nu_res_abs_diff=Cos_Theta_nu_pred-Cos_Theta_nu_true
#Mass squared
Mass_squared_true=((E_nu_true ** 2) - (Px_nu_true ** 2) - (Py_nu_true ** 2) - (Pz_nu_true ** 2)  ) ** 1.0
Mass_squared_pred=((E_nu_pred ** 2) - (Px_nu_pred ** 2) - (Py_nu_pred ** 2) - (Pz_nu_pred ** 2)  ) ** 1.0

#Start plotting

#1D resolutions as a function of true variables
#Energy
binstat.plot_xstat(E_nu_true,E_nu_res_percent,name=plotpath +"/E_nu_res_1D_percent",title="Atmospheric Neutrino Energy Reconstruction",scale='linear',xlabel="True Neutrino Energy",ylabel="Resolution on Neutrino Energy, Percent",range=(0.0,1.01))
plt.close()
binstat.plot_xstat(E_nu_true,E_nu_res_percent_diff,name=plotpath +"/E_nu_res_1D_percent_diff",title="Atmospheric Neutrino Energy Reconstruction",scale='linear',xlabel="True Neutrino Energy",ylabel="Resolution on Neutrino Energy, Percent Difference",range=(0.0,1.01))
plt.close()
binstat.plot_xstat(E_nu_true,E_nu_res_abs_diff,name=plotpath +"/E_nu_res_1D_abs_diff",title="Atmospheric Neutrino Energy Reconstruction",scale='linear',xlabel="True Neutrino Energy",ylabel="Resolution on Neutrino Energy, Absolute Difference",range=(0.0,1.01))
plt.close()
#Cos_Theta
binstat.plot_xstat(Cos_Theta_nu_true,Cos_Theta_nu_res_percent,name=plotpath +"/Cos_Theta_nu_res_1D_percent",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Atmospheric Neutrinos' Cosine of Incoming Angle",xlabel="True Cosine of Incoming Angle",ylabel="Resolution on Cosine of Incoming Angle, Percent")
plt.close()
binstat.plot_xstat(Cos_Theta_nu_true,Cos_Theta_nu_res_percent_diff,name=plotpath +"/Cos_Theta_nu_res_1D_percent_diff",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Atmospheric Neutrinos' Cosine of Incoming Angle",xlabel="True Cosine of Incoming Angle",ylabel="Resolution on Cosine of Incoming Angle, Percent Difference")
plt.close()
binstat.plot_xstat(Cos_Theta_nu_true,Cos_Theta_nu_res_abs_diff,name=plotpath +"/Cos_Theta_nu_res_1D_abs_diff",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Atmospheric Neutrinos' Cosine of Incoming Angle",xlabel="True Cosine of Incoming Angle",ylabel="Resolution on Cosine of Incoming Angle, Absolute Difference")
plt.close()
#Theta
binstat.plot_xstat(Theta_nu_true,Theta_nu_res_percent,name=plotpath +"/Theta_nu_res_1D_percent",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Atmospheric Neutrinos' Incoming Angle",xlabel="True Incoming Angle",ylabel="Resolution on Incoming Angle, Percent")
plt.close()
binstat.plot_xstat(Theta_nu_true,Theta_nu_res_percent_diff,name=plotpath +"/Theta_nu_res_1D_percent_diff",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Atmospheric Neutrinos' Incoming Angle",xlabel="True Incoming Angle",ylabel="Resolution on Incoming Angle, Percent Difference")
plt.close()
binstat.plot_xstat(Theta_nu_true,Theta_nu_res_abs_diff,name=plotpath +"/Theta_nu_res_1D_abs_diff",xrange=(-1.0,1.0),yrange=(-1.0,1.0),title="Atmospheric Neutrinos' Incoming Angle",xlabel="True Incoming Angle",ylabel="Resolution on Cosine of Incoming Angle, Absolute Difference")
plt.close()

#1D resolutions
#Energy
binstat.plot_y_hist(E_nu_res_percent,name=plotpath +"/E_nu_res_percent",range=(-100.0,100.0),bins=500)
plt.close()
binstat.plot_y_hist(E_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_diff",range=(-100.0,100.0),bins=500)
plt.close()
binstat.plot_y_hist(E_nu_res_abs_diff,name=plotpath+"/E_nu_res_abs_diff",range=(-5.0,5.0),bins=500)
plt.close()
#Cos_Theta
binstat.plot_y_hist(Cos_Theta_nu_res_percent,name=plotpath +"/Cos_Theta_nu_res_percent",range=(-100.0,100.0),bins=500)
plt.close()
binstat.plot_y_hist(Cos_Theta_nu_res_percent_diff,name=plotpath +"/Cos_Theta_nu_res_percent_diff",range=(-100.0,100.0),bins=500)
plt.close()
binstat.plot_y_hist(Cos_Theta_nu_res_abs_diff,name=plotpath+"/Cos_Theta_nu_res_abs_diff",range=(-1.0,1.0),bins=500)
plt.close()
#Theta
binstat.plot_y_hist(Theta_nu_res_percent,name=plotpath +"/Theta_nu_res_percent",range=(-100.0,100.0),bins=500)
plt.close()
binstat.plot_y_hist(Theta_nu_res_percent_diff,name=plotpath +"/Theta_nu_res_percent_diff",range=(-100.0,100.0),bins=500)
plt.close()
binstat.plot_y_hist(Theta_nu_res_abs_diff,name=plotpath+"/Theta_nu_res_abs_diff",range=(-180.0,180.0),bins=500)
plt.close()

#1D true and predicted variable distributions for comparisons
#Energy
binstat.plot_y_hist(E_nu_true,E_nu_pred,name=plotpath +"/E_nu_true_vs_pred_1D",range=(-0.01,1.01),bins=500,log=True,labels=["True Energy","Pred. Energy"])
plt.close()
#Cos_Theta
binstat.plot_y_hist(Cos_Theta_nu_true,Cos_Theta_nu_pred,name=plotpath +"/Cos_Theta_nu_true_vs_pred_1D",range=(-1.0,1.0),bins=500,log=True,labels=["True CosTheta","Pred. CosTheta"])
plt.close()
#Theta
binstat.plot_y_hist(Theta_nu_true,Theta_nu_pred,name=plotpath +"/Theta_nu_true_vs_pred_1D",range=(0.,180.0),bins=500,log=True,labels=["True Theta","Pred. Theta"])
plt.close()
#Mass squared
binstat.plot_y_hist(Mass_squared_true,Mass_squared_pred,name=plotpath +"/Mass_squared_nu_true_vs_pred_1D",range=(-0.5,2.0),labels=["True Mass Squared","Pred. Mass Squared"])
plt.close()

#2D plots
#True vs. predicted variables
#Energy
binstat.plot_2d_hist_count(E_nu_true,E_nu_pred,name=plotpath +"/E_nu_2D_true_vs_pred",xrange=(0.0,1.0),yrange=(0.0,1.0),xbins=200,ybins=200,title="Atmospheric Neutrino Energy Reconstruction",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Predicted Neutrino Energy (GeV)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_true,E_nu_pred,name=plotpath +"/E_nu_2D_true_vs_pred_wContours",xrange=(0.0,1.0),yrange=(0.0,1.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Atmospheric Neutrino Energy Reconstruction",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Predicted Neutrino Energy (GeV)",contours=True,contour_labels=True,histogram=True,diag_line=True)
plt.close()

#Cos_Theta
binstat.plot_2d_hist_count(Cos_Theta_nu_true,Cos_Theta_nu_pred,name=plotpath +"/Cos_Theta_nu_2D_true_vs_pred",xrange=(-1.0,1.0),yrange=(-1.0,1.0),xbins=200,ybins=200,title="Atmospheric Neutrinos' Cosine of Incoming Angle (unitless)",scale='linear',zscale='log',xlabel="True Cosine of Incoming Angle (rad)",ylabel="Predicted Cosine of Incoming Angle (unitless)")
plt.close()
binstat.plot_2d_hist_contour(Cos_Theta_nu_true,Cos_Theta_nu_pred,name=plotpath +"/Cos_Theta_nu_2D_true_vs_pred_wContours",xrange=(-1.0,1.0),yrange=(-1.0,1.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Atmospheric Neutrinos' Cosine of Incoming Angle (unitless)",scale='linear',zscale='log',xlabel="True Cosine of Incoming Angle (rad)",ylabel="Predicted Cosine of Incoming Angle (unitless)",contours=True,contour_labels=True,histogram=True,diag_line=True)
plt.close()

#Theta
binstat.plot_2d_hist_count(Theta_nu_true,Theta_nu_pred,name=plotpath +"/Theta_nu_2D_true_vs_pred",xrange=(0.,180.0),yrange=(0.,180.0),xbins=200,ybins=200,title="Atmospheric Neutrinos' Incoming Angle",scale='linear',zscale='log',xlabel="True Incoming Angle (deg)",ylabel="Predicted Incoming Angle (deg)")
plt.close()
binstat.plot_2d_hist_contour(Theta_nu_true,Theta_nu_pred,name=plotpath +"/Theta_nu_2D_true_vs_pred_wContours",xrange=(0.,180.0),yrange=(0.,180.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Atmospheric Neutrinos' Incoming Angle",scale='linear',zscale='log',xlabel="True Incoming Angle (deg)",ylabel="Predicted Incoming Angle (deg)",contours=True,contour_labels=True,histogram=True,diag_line=True)
plt.close()

#Energy resolutions as a function of Theta resolutions, all different types
binstat.plot_2d_hist_count(E_nu_res_percent,Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_vs_Theta_nu_res_percent",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="Theta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent,Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_vs_Theta_nu_res_percent_wContours",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="Theta Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True,circle=True,radius=10)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_vs_Theta_nu_res_percent_diff",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="Theta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_vs_Theta_nu_res_percent_diff_wContours",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="Theta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True,circle=True,radius=10)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_vs_Theta_nu_res_abs_diff",xrange=(-100.0,100.0),yrange=(-180.0,180.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="Theta Resolution, Absolute Difference (Degrees)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_vs_Theta_nu_res_abs_diff_wContours",xrange=(-100.0,100.0),yrange=(-180.0,180.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="Theta Resolution, Absolute Difference (Degrees)",contours=True,contour_labels=True,histogram=True,circle=True,radius=10)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent_diff,Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_diff_vs_Theta_nu_res_percent",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="Theta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent_diff,Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_diff_vs_Theta_nu_res_percent_wContours",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="Theta Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True,circle=True,radius=10)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent_diff,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Theta_nu_res_percent_diff",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="Theta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent_diff,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Theta_nu_res_percent_diff_wContour",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="Theta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True,circle=True,radius=10)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent_diff,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Theta_nu_res_abs_diff",xrange=(-100.0,100.0),yrange=(-180.0,180.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="Theta Resolution, Absolute Difference (Degrees)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent_diff,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Theta_nu_res_abs_diff_wContour",xrange=(-100.0,100.0),yrange=(-180.0,180.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="Theta Resolution, Absolute Difference (Degrees)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10,y_semiminor=30)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_abs_diff,Theta_nu_res_percent,name=plotpath +"/E_nu_res_abs_diff_vs_Theta_nu_res_percent",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="Theta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_abs_diff,Theta_nu_res_percent,name=plotpath +"/E_nu_res_abs_diff_vs_Theta_nu_res_percent_wContour",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="Theta Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=0.250,y_semiminor=10)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_abs_diff,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Theta_nu_res_percent_diff",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="Theta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_abs_diff,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Theta_nu_res_percent_diff_wContour",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="Theta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=0.250,y_semiminor=10)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_abs_diff,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Theta_nu_res_abs_diff",xrange=(-1.2,1.2),yrange=(-180.0,180.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="Theta Resolution, Absolute Difference (Degrees)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_abs_diff,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Theta_nu_res_abs_diff_wContour",xrange=(-1.2,1.2),yrange=(-180.0,180.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="Theta Resolution, Absolute Difference (Degrees)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=0.250,y_semiminor=30)
plt.close()



#Energy resolutions as a function of Cos_Theta resolutions, all different types
binstat.plot_2d_hist_count(E_nu_res_percent,Cos_Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_vs_Cos_Theta_nu_res_percent",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="CosTheta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent,Cos_Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_vs_Cos_Theta_nu_res_percent_wContour",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="CosTheta Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=20.0)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent,Cos_Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_vs_Cos_Theta_nu_res_percent_diff",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="CosTheta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent,Cos_Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_vs_Cos_Theta_nu_res_percent_diff_wContour",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="CosTheta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=20.0)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent,Cos_Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_vs_Cos_Theta_nu_res_abs_diff",xrange=(-100.0,100.0),yrange=(-1.0,1.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="CosTheta Resolution, Absolute Difference")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent,Cos_Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_vs_Cos_Theta_nu_res_abs_diff_wContour",xrange=(-100.0,100.0),yrange=(-1.0,1.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent (%)",ylabel="CosTheta Resolution, Absolute Difference",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=0.2)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent_diff,Cos_Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_diff_vs_Cos_Theta_nu_res_percent",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="CosTheta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent_diff,Cos_Theta_nu_res_percent,name=plotpath +"/E_nu_res_percent_diff_vs_Cos_Theta_nu_res_percent_wContour",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="CosTheta Resolution, Percent( %)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=20.0)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent_diff,Cos_Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Cos_Theta_nu_res_percent_diff",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="CosTheta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent_diff,Cos_Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Cos_Theta_nu_res_percent_diff_wContour",xrange=(-100.0,100.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="CosTheta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=20.0)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_percent_diff,Cos_Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Cos_Theta_nu_res_abs_diff",xrange=(-100.0,100.0),yrange=(-2.0,2.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="CosTheta Resolution, Absolute Difference")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_percent_diff,Cos_Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_percent_diff_vs_Cos_Theta_nu_res_abs_diff_wContour",xrange=(-100.0,100.0),yrange=(-2.0,2.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Percent Difference (%)",ylabel="CosTheta Resolution, Absolute Difference",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=0.2)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_abs_diff,Cos_Theta_nu_res_percent,name=plotpath +"/E_nu_res_abs_diff_vs_Cos_Theta_nu_res_percent",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="CosTheta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_abs_diff,Cos_Theta_nu_res_percent,name=plotpath +"/E_nu_res_abs_diff_vs_Cos_Theta_nu_res_percent_wContour",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="CosTheta Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=0.25,y_semiminor=20.0)
plt.close()
                           
binstat.plot_2d_hist_count(E_nu_res_abs_diff,Cos_Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Cos_Theta_nu_res_percent_diff",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="CosTheta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_abs_diff,Cos_Theta_nu_res_percent_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Cos_Theta_nu_res_percent_diff_wContour",xrange=(-1.2,1.2),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="CosTheta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=0.25,y_semiminor=20.0)
plt.close()

binstat.plot_2d_hist_count(E_nu_res_abs_diff,Cos_Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Cos_Theta_nu_res_abs_diff",xrange=(-1.2,1.2),yrange=(-2.0,2.0),xbins=200,ybins=200,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="CosTheta Resolution, Absolute Difference")
plt.close()
binstat.plot_2d_hist_contour(E_nu_res_abs_diff,Cos_Theta_nu_res_abs_diff,name=plotpath +"/E_nu_res_abs_diff_vs_Cos_Theta_nu_res_abs_diff",xrange=(-1.2,1.2),yrange=(-2.0,2.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Overall Performance of Kinematic Reconstruction",scale='linear',zscale='log',xlabel="Energy Resolution, Absolute Difference (GeV)",ylabel="CosTheta Resolution, Absolute Difference",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=0.25,y_semiminor=0.2)
plt.close()




#Energy resolution vs true energy
binstat.plot_2d_hist_count(E_nu_true,E_nu_res_percent,name=plotpath +"/E_nu_true_vs_E_nu_res_percent",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Energy Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Energy Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_true,E_nu_res_percent,name=plotpath +"/E_nu_true_vs_E_nu_res_percent_wContour",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Energy Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Energy Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True)
plt.close()

binstat.plot_2d_hist_count(E_nu_true,E_nu_res_percent_diff,name=plotpath +"/E_nu_true_vs_E_nu_res_percent_diff",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Energy Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Energy Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_true,E_nu_res_percent_diff,name=plotpath +"/E_nu_true_vs_E_nu_res_percent_diff_wContour",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Energy Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Energy Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True)
plt.close()

binstat.plot_2d_hist_count(E_nu_true,E_nu_res_abs_diff,name=plotpath +"/E_nu_true_vs_E_nu_res_abs_diff",xrange=(0.0,1.0),yrange=(-1.0,1.0),xbins=200,ybins=200,title="Energy Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Energy Resolution, Absolute Difference")
plt.close()
binstat.plot_2d_hist_contour(E_nu_true,E_nu_res_abs_diff,name=plotpath +"/E_nu_true_vs_E_nu_res_abs_diff_wContour",xrange=(0.0,1.0),yrange=(-1.0,1.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Energy Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Energy Resolution, Absolute Difference",contours=True,contour_labels=True,histogram=True)
plt.close()






#Theta resolution vs true energy
binstat.plot_2d_hist_count(E_nu_true,Theta_nu_res_percent,name=plotpath +"/E_nu_true_vs_Theta_nu_res_percent",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Theta Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Theta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_true,Theta_nu_res_percent,name=plotpath +"/E_nu_true_vs_Theta_nu_res_percent_wContour",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Theta Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Theta Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True)
plt.close()

binstat.plot_2d_hist_count(E_nu_true,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_true_vs_Theta_nu_res_percent_diff",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Theta Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Theta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(E_nu_true,Theta_nu_res_percent_diff,name=plotpath +"/E_nu_true_vs_Theta_nu_res_percent_diff_wContour",xrange=(0.0,1.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Theta Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Theta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True)
plt.close()

binstat.plot_2d_hist_count(E_nu_true,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_true_vs_Theta_nu_res_abs_diff",xrange=(0.0,1.0),yrange=(-180.0,180.0),xbins=200,ybins=200,title="Theta Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Theta Resolution, Absolute Difference")
plt.close()
binstat.plot_2d_hist_contour(E_nu_true,Theta_nu_res_abs_diff,name=plotpath +"/E_nu_true_vs_Theta_nu_res_abs_diff_wContour",xrange=(0.0,1.0),yrange=(-180.0,180.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Theta Reconstruction Performance Across Energies",scale='linear',zscale='log',xlabel="True Neutrino Energy (GeV)",ylabel="Theta Resolution, Absolute Difference",contours=True,contour_labels=True,histogram=True)
plt.close()







#Theta resolution vs true theta
binstat.plot_2d_hist_count(Theta_nu_true,Theta_nu_res_percent,name=plotpath +"/Theta_nu_true_vs_Theta_nu_res_percent",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Theta Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Theta Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(Theta_nu_true,Theta_nu_res_percent,name=plotpath +"/Theta_nu_true_vs_Theta_nu_res_percent_wContour",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Theta Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Theta Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=20.0)
plt.close()

binstat.plot_2d_hist_count(Theta_nu_true,Theta_nu_res_percent_diff,name=plotpath +"/Theta_nu_true_vs_Theta_nu_res_percent_diff",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Theta Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Theta Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(Theta_nu_true,Theta_nu_res_percent_diff,name=plotpath +"/Theta_nu_true_vs_Theta_nu_res_percent_diff_wContour",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Theta Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Theta Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True,ellipse=True,x_semimajor=10.0,y_semiminor=20.0)
plt.close()

binstat.plot_2d_hist_count(Theta_nu_true,Theta_nu_res_abs_diff,name=plotpath +"/Theta_nu_true_vs_Theta_nu_res_abs_diff",xrange=(0.0,180.0),yrange=(-180.0,180.0),xbins=200,ybins=200,title="Theta Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Theta Resolution, Absolute Difference")
plt.close()
binstat.plot_2d_hist_contour(Theta_nu_true,Theta_nu_res_abs_diff,name=plotpath +"/Theta_nu_true_vs_Theta_nu_res_abs_diff_wContour",xrange=(0.0,180.0),yrange=(-180.0,180.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Theta Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Theta Resolution, Absolute Difference",contours=True,contour_labels=True,histogram=True)
plt.close()






#Energy resolution vs true theta
binstat.plot_2d_hist_count(Theta_nu_true,E_nu_res_percent,name=plotpath +"/Theta_nu_true_vs_E_nu_res_percent",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Energy Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Energy Resolution, Percent (%)")
plt.close()
binstat.plot_2d_hist_contour(Theta_nu_true,E_nu_res_percent,name=plotpath +"/Theta_nu_true_vs_E_nu_res_percent_wContour",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Energy Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Energy Resolution, Percent (%)",contours=True,contour_labels=True,histogram=True)
plt.close()

binstat.plot_2d_hist_count(Theta_nu_true,E_nu_res_percent_diff,name=plotpath +"/Theta_nu_true_vs_E_nu_res_percent_diff",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,title="Energy Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Energy Resolution, Percent Difference (%)")
plt.close()
binstat.plot_2d_hist_contour(Theta_nu_true,E_nu_res_percent_diff,name=plotpath +"/Theta_nu_true_vs_E_nu_res_percent_diff_wContour",xrange=(0.0,180.0),yrange=(-100.0,100.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Energy Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Energy Resolution, Percent Difference (%)",contours=True,contour_labels=True,histogram=True)
plt.close()

binstat.plot_2d_hist_count(Theta_nu_true,E_nu_res_abs_diff,name=plotpath +"/Theta_nu_true_vs_E_nu_res_abs_diff",xrange=(0.0,180.0),yrange=(-180.0,180.0),xbins=200,ybins=200,title="Energy Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Energy Resolution, Absolute Difference")
plt.close()
binstat.plot_2d_hist_contour(Theta_nu_true,E_nu_res_abs_diff,name=plotpath +"/Theta_nu_true_vs_E_nu_res_abs_diff_wContour",xrange=(0.0,180.0),yrange=(-180.0,180.0),xbins=200,ybins=200,xbins_contour=25,ybins_contour=25,title="Energy Reconstruction Performance Across Angles",scale='linear',zscale='log',xlabel="True Neutrino Theta (Degree)",ylabel="Energy Resolution, Absolute Difference",contours=True,contour_labels=True,histogram=True)
plt.close()



