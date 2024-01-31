# ****************************************************************************
# Name       : utils.py
# Author     : Andres Valdez
# Version    : 1.0
# Description: User defined tools for data analysis
# Data	     : 01-21-2024
# ****************************************************************************

from __future__ import unicode_literals
import os, sys
import numpy as np
import scipy as scp
import pandas as pd
from scipy import optimize
from scipy.signal import butter,filtfilt
import time

from statsmodels.tsa.seasonal import seasonal_decompose

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
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    window_size = np.abs(int(window_size))
    order = np.abs(int(order))
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def butter_lowpass_filter(data, cutoff, fs, order):
    """
    source https://medium.com/analytics-vidhya/how-to-filter-noise-with-a-low-pass-filter-python-885223e5e9b7
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def read_csv_data(filename='20230920_001_DarkFullData.csv'):
    # Read the data
    data = np.loadtxt(filename, delimiter=',', usecols=(0, 1), unpack=True,skiprows=2)
    print('Analyzed file  :',filename)
    
    # time X Opt.Data
    t_data = data[0]
    dp_val = data[1]
    
    #print(t_data.shape,dp_val.shape)
    
    return t_data, dp_val
    
def get_spectrum(t,y,label,kval=10,full_output=False):
    """
    This function will return the first kval vibration modes of the glottal flow signal
    t         :: Time array
    y         :: Signal to analyze
    kval      :: Number of modes adopted
    """
    y_fft  = scp.fftpack.fft(y)
    y_amp  = 2.0 / t.size * np.abs(y_fft)
    y_freq = np.abs(scp.fftpack.fftfreq(t.size))
    
    # Get the amplitudes
    sa     = pd.Series(y_amp).nlargest(kval).round(3).astype(float).tolist()
    
    # Get the frequencies
    magnitudes = abs(y_fft[np.where(y_freq >= 0)])
    sf         = np.sort((np.argpartition(magnitudes, -kval)[-kval:])/t[-1])
    
    print('#############################################################')
    print('Signal:',label)
    print('Main freqs.',sf)
    print('Main amps. ',sa)
    print('#############################################################')
    print('')
    
    if(full_output):
        return y_freq , y_amp, np.array(sa), np.array(sf)
    else:
        return y_freq, y_amp

def get_derivative(t,y):
    """
    This function returns the time derivative  of y (Re-shape-ing the arrays)
    """
    
    dy = np.gradient(y)
    dt = np.gradient(t)
    
    return dy / dt

def get_shocklocs(t,dp):
    """
    This function returns the shock locations for the pressure drop
    """
    # Get the derivative of pressure drop quasi-convoluted with dp
    dertest   = get_derivative(t,dp) * dp / (np.linalg.norm(dp))
    
    # Filter the derivative
    dertest   =  savitzky_golay(dertest, 70, 2)
    
    # This is a trick, to evaluate regions of common-jumps
    index_foo = dp - dertest
    
    longidx = np.argwhere(np.diff(np.sign(dp - index_foo))).flatten()
    
    # Now remove the trivial cases
    idx = []
    pval = 0.01
    for k in longidx:
        if( dp[k] > pval):
            idx.append(k)
            pval = dp[k]
    
    return idx

########################################################################
# Here comes the main function
########################################################################

if __name__ == "__main__":

    t0 = time.process_time() # Here start count time

    colors = ['red','green','blue']
    labels = ['data','Butter','Savitzky-Golay']
    
    # Get the data
    t_data , dp_val = read_csv_data(filename='20230920_001_ShockData.csv')
    
    t_data , dp_val = t_data[:15000] , dp_val[:15000]
    
    # Filter the data
    dp_savgol =  savitzky_golay(dp_val, 50, 2)
    dp_butter = butter_lowpass_filter(dp_savgol,2,300,2)
    
    # Now I have tow stagered filters, but will improve merging them
    test      = np.sqrt(np.abs(dp_butter*dp_savgol))
    
    #period = int(0.5*t_data[-1] / (t_data[2]-t_data[1]))
    #results = seasonal_decompose(dp_val, model='additive', period=period)
    #test = results.trend

    # Get vibration pattern
    fx , ampx = get_spectrum(t_data , dp_val, labels[0])
    fx_svg , ampx_svg = get_spectrum(t_data , dp_savgol,labels[1])
    fx_btt , ampx_btt = get_spectrum(t_data , dp_butter,labels[2])
    fx_tst , ampx_tst = get_spectrum(t_data , test,'test')
    
        
    fig , axs = plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    
    axs[0].plot(t_data,dp_val,'-',lw=2,color=colors[0],label=labels[0])
    #axs[0].plot(t_data,dp_butter,'-',lw=2,color=colors[1],label=labels[1])
    #axs[0].plot(t_data,dp_savgol,'-',lw=2,color=colors[2],label=labels[2])

    axs[0].plot(t_data,test,'-',lw=2,color=colors[1],label='test')


    axs[1].semilogy(fx,ampx,'-',lw=2,color=colors[0])
    #axs[1].semilogy(fx_svg,ampx_svg,'-',lw=2,color=colors[1])
    #axs[1].semilogy(fx_btt,ampx_btt,'-',lw=2,color=colors[2])

    axs[1].semilogy(fx_tst,ampx_tst,'-',lw=2,color=colors[1])

    axs[0].set_xlabel('time (ms)')
    axs[0].set_ylabel('$\delta_p$')

    axs[1].set_xlabel('Freq.')
    axs[1].set_ylabel('amplitude')    
    
    axs[0].legend(loc='best')
    
    plt.show()
    
    t1 = time.process_time() # Here end counting time
    
    print("Elapsed time to solve: %0.2f secs." % (t1-t0))



