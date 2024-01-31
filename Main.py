# ****************************************************************************
# Name       : Main.py
# Author     : Andres Valdez
# Version    : 1.0
# Description: Rapid test function
# Data	     : 01-21-2024
# ****************************************************************************

from utils import *

########################################################################
# Here comes the main testing function
########################################################################

if __name__ == "__main__":

    t0 = time.process_time() # Here start count time

    colors = ['red','green','blue']
    labels = ['data','Butter','Savitzky-Golay']
    
    # Get the data
    t_data , dp_val = read_csv_data(filename='20230920_001_ShockData.csv')
    
    t_data , dp_val = t_data[8500:12000] , dp_val[8500:12000]
    
    # Filter the data
    dp_savgol =  savitzky_golay(dp_val, 50, 2)
    dp_butter = butter_lowpass_filter(dp_savgol,2,300,2)
    
    # Now I have tow stagered filters, but will improve merging them
    test      = np.sqrt(np.abs(dp_butter*dp_savgol))
    
    idx = get_shocklocs(t_data,test)
    
    # Now I get the derivative
    #dertest   = get_derivative(t_data,test) * test / (np.linalg.norm(test))
    
    #dertest   =  savitzky_golay(dertest, 70, 2)
    
    #index_foo = test - dertest
    
    #longidx = np.argwhere(np.diff(np.sign(test - index_foo))).flatten()
    
    # Now remove the trivial cases
    #idx = []
    #for k in longidx:
    #    if(test[k] > 0.015):
    #        idx.append(k)
    
    fig , axs = plt.subplots(1,1,figsize=(12,6),constrained_layout=True)
    
    axs.plot(t_data,dp_val,'-',lw=2,color=colors[0],label='Raw data')
    
    axs.plot(t_data,test,'-',lw=2,color=colors[1],label='Filtered')
    
    axs.plot(t_data[idx], test[idx], 'o', color=colors[2], label='potential shocks')
    
    #axs.plot(t_data,dertest,'-',lw=2,color=colors[1],label='der test (scaled)')

    #axs.plot(t_data,test - dertest,'-',lw=2,color=colors[2],label='index')

    #axs.plot(t_data,(test - dertest),'-',lw=2,color=colors[2],label='index (scaled)')
    
    axs.set_xlabel('time (ms)')
    axs.set_ylabel('$\delta_p$')
    
    axs.legend(loc='best')
    
    #plt.savefig('TestRes.png')
    
    plt.show()
    
    sys.exit()
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



