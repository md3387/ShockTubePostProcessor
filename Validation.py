# ****************************************************************************
# Name       : optical_models.py
# Author     : Andres Valdez
# Version    : 1.0
# Description: Several mathematical models for bio screen analysis
# Data	     : 10-11-2021
# ****************************************************************************

from __future__ import unicode_literals
import os, sys
import numpy as np
import scipy as scp
import pandas as pd
from scipy import optimize
from scipy.signal import butter,filtfilt
import time

from Main import *

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

#For 0.45\textwidth figs works ok
mpl.rcParams['axes.labelsize']  = 17
mpl.rcParams['xtick.labelsize'] = 17
mpl.rcParams['ytick.labelsize'] = 17
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['text.usetex']     = True
mpl.rcParams['font.family']     = 'serif'

import warnings
warnings.filterwarnings("ignore")

np.random.seed(123)

########################################################################
# Routines and functions
########################################################################
    
def read_xlsx_data(filename='validation.xlsx',kval=1):
    """
    This function will upload the content of the excel file to validate the filtering
    """
    
    data_pd = pd.read_excel(filename, skiprows=4)
    t_data  = data_pd.iloc[:,0].to_numpy()
    dp_val  = data_pd.iloc[:,kval].to_numpy()
     
    return t_data , dp_val

########################################################################
# Here comes the main function
########################################################################

if __name__ == "__main__":

    t0 = time.process_time() # Here start count time

    colors = ['red','green','blue']
    labels = ['data','Butter','Savitzky-Golay']
    
    init = [0.8,1.9,1.9,1.9,1.9,1.9,1.9,1.9,1.9,0.8,1.9,1.9,1.9]
    end  = [2.0,2.5,2.5,2.5,2.5,2.5,2.5,2.7,2.5,1.2,2.5,2.7,3.0]
    
    for k in range(13):
        # Get the data
        t_data , dp_val = read_xlsx_data('validation.xlsx',k+1)
        
        #t_data , dp_val = t_data[:15000] , dp_val[:15000]
    
        # Filter the data
        dp_savgol =  savitzky_golay(dp_val, 50, 2)
        dp_butter = butter_lowpass_filter(dp_savgol,2,300,2)
    
        # Now I have tow stagered filters, but will improve merging them
        test      = np.sqrt(np.abs(dp_butter*dp_savgol))
    
        #period = int(0.5*t_data[-1] / (t_data[2]-t_data[1]))
        #results = seasonal_decompose(dp_val, model='additive', period=period)
        #test = results.trend

        fig , axs = plt.subplots(1,1,figsize=(12,10),constrained_layout=True)
    
        axs.plot(t_data,dp_val,'-',lw=2,color=colors[0],label=labels[0])
        #axs[0].plot(t_data,dp_butter,'-',lw=2,color=colors[1],label=labels[1])
        #axs[0].plot(t_data,dp_savgol,'-',lw=2,color=colors[2],label=labels[2])

        axs.plot(t_data,test,'-',lw=2,color=colors[1],label='test')


        axs.set_xlim([init[k],end[k]])
        axs.set_xlabel('time (ms)')
        axs.set_ylabel('$\delta_p$')
        axs.legend(loc='best')
    
        plt.savefig('Filtered/Solution_' + str(k+1) + '.png')

        print('Done with Column:',k+1)
    
    t1 = time.process_time() # Here end counting time
    
    print("Elapsed time to solve: %0.2f mins." % ((t1-t0)/60.0))



