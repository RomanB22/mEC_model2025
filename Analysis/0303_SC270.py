## AUM SRI SAI RAM


import numpy as np
from scipy import signal
from scipy.signal import cwt, morlet2, welch
import matplotlib.pyplot as plt
from funcs import *
import glob, os
from netpyne.analysis.tools import loadData
from itertools import compress
import matplotlib.gridspec as gridspec

##############################################
##### LOAD THE NETPYNE PICKLE FILE############

cellTrace={}
VC = 0 #Clamped voltage
rs=1e-4
Factor = 1e3
for file in glob.glob("../0227_Pastoll/Ant&Ca_0_Kv3_1.0_Kv7_1.0_gsinInh7.0_gsinExc0.0_*_II_0.0_GAP_F_*.pkl"):
    print('Loading the file from specified folder')
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][13:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    cellTrace['0'] = fileInfo['simData']['V_soma']['cell_270']#20
    tSim = fileInfo['simData']['t']

cellTr0 = np.array(cellTrace['0'])
FileID = "0306_Pastoll/SC270_" 
############################################################################################
"""Here is the part to modifiy according to the simulation configuration/experimental setup"""
dt = 1e-2 # in ms
CycToPlot = 2 # Plot the last CycToPlot cycles
fTheta = 8
thetaPeriod = 1000/fTheta
maxPower= 10000 #2 # Sets the maximum for the scalogram. 2 for f120 signal, 7 for f95 and notch

############################################################################################

CapCurrentLeft = 70; CapCurrentRight = 80;
waveletPlotPoints = 0
y0 = (VC-cellTr0)/rs*Factor
y0 -= y0[-1] #offsetting based on zero
try:
    index = int(np.argwhere(y0<minV)[0])
except:
    index = CapCurrentLeft
print('Printing legths of cellTr, tSim, and y0 : ',len(cellTr0),len(tSim),len(y0))

b, a = signal.butter(3, [50, 300], fs=1/(dt*1e-3), btype='band')
rates = signal.filtfilt(b, a, y0) #Here is where we are using the band pass filter and using the signal to get gamma waveform
rates -= rates[-1] #again, offsetting based on zero
print('Length of the assigned array is',np.shape(rates))

############################################################################################

dt_rec = dt; # Down sample the rates (in ms)

fs = np.arange(50.,300.1,0.2);

chunkLength = int(np.round(thetaPeriod/dt_rec));

CycToPlot = 8


yav=np.zeros((1,chunkLength))
for j in range(CycToPlot):
    yav += rates[(CycToPlot-j)*chunkLength-waveletPlotPoints:(CycToPlot-j+1)*chunkLength-waveletPlotPoints]#We are plotting the average here, similar to Brandon's
yav = yav[0]/CycToPlot

print('Printing the length of rates and average and checking if it is correct: ', len(yav))

print('Trying all this stuff to generate average wavelets, Lets see!!!!!')

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
Power_Index = 2809010.088019215
plotPow(cwtPowAvg,Ncycs=1,dt=dt_rec,fs=fs,levels=np.linspace(0.,Power_Index,150,endpoint=True))
plt.ylabel("Gamma Frequency (Hz)")
plt.xlabel("Stim Theta Phase (rad)")


plt.savefig(FileID+"Pow_AVG.eps", bbox_inches='tight')
plt.savefig(FileID+"Pow_AVG.png", bbox_inches='tight')
plt.savefig(FileID+"Pow_AVG.svg", bbox_inches='tight')
plt.savefig(FileID+"Pow_AVG.pdf", bbox_inches='tight')
#plt.show()
plt.clf();


############################################################################################
#### Trying All Together at Once based on Older Figures - Let me see


fig = plt.figure(figsize=(10, 8))
outer = gridspec.GridSpec(1, 1, wspace=0.3, hspace=0.1)

inner = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.1, height_ratios=[1,1,1])

ax0 = plt.Subplot(fig, inner[0])
ax0.plot(tSim, y0, 'b-',label='Unfiltered')
ax0.set_xlim(sim_time-4.*125.,sim_time);
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(True)
ax0.set_ylabel("Current (pA)")
#ax0.set_ylim(1.2*min(y0[-1000:-1]),1.2*max(y0[-1000:-1]))
ax0.set_ylim([-100,1000])
ax0.set_xticks([])
ax0.set_xticklabels([])
ax0.legend(loc='upper right')
fig.add_subplot(ax0)

ax1 = plt.Subplot(fig, inner[1])
ax1.plot(tSim, rates, 'r-',label='Filtered 50-300 hz')
ax1.set_xlim(sim_time-4.*125.,sim_time);
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(True)
ax1.set_ylabel("Current (pA)")
#ax1.set_ylim(1.2*min(y1[-1000:-1]),1.2*max(y1[-1000:-1]))
ax1.set_ylim([-450,450])
ax1.set_xticks([])
ax1.set_xticklabels([])
ax1.legend(loc='upper right')
fig.add_subplot(ax1)

gsin = 5

ax2 = plt.Subplot(fig, inner[2])
t = np.linspace(sim_time-4.*125.,sim_time,4*125+1)
ax2.set_xlim(sim_time-4.*125.,sim_time)
ax2.plot(t, gsin*np.sin(2*np.pi*fTheta*t*1e-3-np.pi/2),color = 'k')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.set_xticks([])
ax2.set_yticks([])
fig.add_subplot(ax2)



fig.tight_layout()
fig.savefig(FileID+'Currents.png', bbox_inches='tight')
fig.savefig(FileID+'Currents.eps', bbox_inches='tight')

