## AUM SRI SAI RAM


import numpy as np
from scipy import signal
from scipy.signal import cwt, morlet2, welch
import matplotlib.pyplot as plt
from funcs_theta import *
import glob, os
from netpyne.analysis.tools import loadData
from itertools import compress
import matplotlib.gridspec as gridspec

FileID = "0306_Plots/Testing_Sine_" 
############################################################################################
"""Here is the part to modifiy according to the simulation configuration/experimental setup"""
dt = 1e-2 # in ms
CycToPlot = 2 # Plot the last CycToPlot cycles
fTheta = 8
thetaPeriod = 1000/fTheta
maxPower= 10000 #2 # Sets the maximum for the scalogram. 2 for f120 signal, 7 for f95 and notch
waveletPlotPoints = 0
############################################################################################
sim_time = 12*125
t = np.linspace(sim_time-125,sim_time,12*12500+1)
gsin = 1
sine_wave = gsin*np.sin(2*np.pi*fTheta*t*1e-3-np.pi/2)
print('Printing Sim_Time and length of sine_wave: ',sim_time,len(sine_wave))

############################################################################################

dt_rec = dt; # Down sample the rates (in ms)

fs = np.arange(0.,20.1,0.1);

chunkLength = int(np.round(thetaPeriod/dt_rec));

print('Printing length of fs and value of chunkLength: ', len(fs),chunkLength)

CycToPlot = 1


yav=np.zeros((1,chunkLength))
for j in range(CycToPlot):
    yav += sine_wave[(CycToPlot-j)*chunkLength:(CycToPlot-j+1)*chunkLength-waveletPlotPoints]#We are plotting the average here, similar to Brandon's
yav = yav[0]/CycToPlot

print('Printing the length of rates and average and checking if it is correct: ', len(yav))

frequencies, psd = calculate_psd(yav, dt_rec)

# Find peak frequency
peak_frequency = find_peak_frequency(frequencies, psd)

# Plot PSD
plt.figure(figsize=(10, 6))
plt.plot(frequencies, psd)
plt.title('Power Spectral Density (PSD)')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power/Frequency (dB/Hz or similar)')
plt.xlim([0,300])
plt.axvline(peak_frequency, color='red', linestyle='--', label=f'Peak: {peak_frequency:.2f} Hz') # Add a line at the peak frequency
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig(FileID+"PSD.eps", dpi=300)
plt.savefig(FileID+"PSD.png", dpi=300)
plt.clf();
# Print peak frequency
print(f"Peak frequency in PSD: {peak_frequency:.2f} Hz")


cwtPowAvg = compPWT(yav,dt=dt_rec*1e-3,fs=fs);
maxPowerAvg = np.max(np.max(cwtPowAvg))
print('Maximum Power in the average case is',maxPowerAvg)

# Find the indices of the maximum power
max_power_index = np.unravel_index(np.argmax(cwtPowAvg), cwtPowAvg.shape)
max_freq_index, max_time_index = max_power_index

# Get the frequency and time at the maximum power
max_freq = fs[max_freq_index]
max_time = max_time_index * dt_rec  # Time in milliseconds

print(f"Maximum power occurs at frequency: {max_freq} Hz")
print(f"Maximum power occurs at time: {max_time} ms (or {max_time * 1e-3} s)")
Power_Index = maxPowerAvg
plotPow(cwtPowAvg,Ncycs=1,dt=dt_rec,fs=fs,levels=np.linspace(0.,Power_Index,150,endpoint=True))
plt.ylabel("Gamma Frequency (Hz)")
plt.xlabel("Stim Theta Phase (rad)")


plt.savefig(FileID+"Pow_AVG.eps", bbox_inches='tight')
plt.savefig(FileID+"Pow_AVG.png", bbox_inches='tight')
plt.savefig(FileID+"Pow_AVG.svg", bbox_inches='tight')
plt.savefig(FileID+"Pow_AVG.pdf", bbox_inches='tight')
#plt.show()
plt.clf();

