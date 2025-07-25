import numpy as np
from funcsAux import *
import matplotlib.pyplot as plt
import glob, os
from netpyne.analysis.tools import loadData
from itertools import compress
import matplotlib.gridspec as gridspec
gsin = 5
f = 8
dt_rec = 1e-2 # Down sample the rates
minFreq=30.
maxFreq=250.
maxPower= 2500 # 100. # Sets the maximum for the scalogram. Varies according to the power of the signal. Change it for better visualization
colorbar=True
numBins = 2*125+1
###################################
fs = np.arange(minFreq,maxFreq,3.)
levelsWav =  np.linspace(0., maxPower, 20, endpoint=True)

def compute_full_autocorrelation(data):
    n = len(data)
    data_mean = np.mean(data)
    centered_data = data - data_mean
    autocorr = np.correlate(centered_data, centered_data, mode='full')
    return autocorr

i = 0
spktFSAux = {}

fileID = "0303_Pastoll/0303_"

for file in glob.glob("../0218_With_IPSCs/Ant*gsinInh7.0*gsinExc0.0*Hyper*_Noise_F*.pkl"):
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][6:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    maskFS = np.array(fileInfo['simData']['spkid'])< 100
    spktFSAux['0'] = np.array( list( compress(fileInfo['simData']['spkt'], maskFS) ) )

spktFS0 = spktFSAux['0']

for file in glob.glob("../0227_Pastoll/Ant&Ca_0_Kv3_1.0_Kv7_1.0_gsinInh7.0_gsinExc0.0_EI_0.3_IE_0.2_II_0.0_GAP_F_VC0.0Hyper_Noise_F*.pkl"):
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][6:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    maskFS = np.array(fileInfo['simData']['spkid'])< 100
    spktFSAux['1'] = np.array( list( compress(fileInfo['simData']['spkt'], maskFS) ) )
spktFS1 = spktFSAux['1']

for file in glob.glob("../0303_Huang/Ant&Ca_0_Kv3_1.0_Kv7_1.0_gsinInh7.0_gsinExc0.0_EI_0.3_IE_0.2_II_0.0_GAP_T_VC0.0Hyper_Noise_F*.pkl"):
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][8:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    maskFS = np.array(fileInfo['simData']['spkid'])< 100
    spktFSAux['2'] = np.array( list( compress(fileInfo['simData']['spkt'], maskFS) ) )
    i += 1
spktFS2 = spktFSAux['2']

histBins = np.linspace(sim_time-2.*125.,sim_time,numBins)
#histFS0 = np.histogram(spktFS0,histBins)
histFS0 = plt.hist(spktFS0, np.linspace(sim_time-2.*125.,sim_time, numBins),color = 'k', label='FS. Stimulate SC Only - Hyper')
plt.xlim(sim_time-2.*125.,sim_time)
plt.ylim(0.,50.)
#plt.hist(histFS0, bins = histBins)
#plt.show()

hist_counts0 = histFS0[0]


histFS1 = plt.hist(spktFS1, np.linspace(sim_time-2.*125.,sim_time, numBins),color = 'k', label='FS. Stimulate SC Only - Hyper')
plt.xlim(sim_time-2.*125.,sim_time)
plt.ylim(0.,50.)
#plt.hist(histFS0, bins = histBins)
#plt.show()

hist_counts1 = histFS1[0]

histFS2 = plt.hist(spktFS2, np.linspace(sim_time-2.*125.,sim_time, numBins),color = 'k', label='FS. Stimulate SC Only - Hyper')
plt.xlim(sim_time-2.*125.,sim_time)
plt.ylim(0.,50.)
#plt.hist(histFS0, bins = histBins)
#plt.show()

hist_counts2 = histFS2[0]

# Compute the full auto-correlation of the histogram counts
full_autocorr_Hyper = compute_full_autocorrelation(hist_counts0)
full_autocorr_Shunt = compute_full_autocorrelation(hist_counts1)
full_autocorr_Uniform = compute_full_autocorrelation(hist_counts2)

# Normalize by the value at lag zero (middle of the array)
mid_point_Hyper = len(full_autocorr_Hyper) // 2
normalized_autocorr_Hyper = full_autocorr_Hyper / full_autocorr_Hyper[mid_point_Hyper]

mid_point_Shunt = len(full_autocorr_Shunt) // 2
normalized_autocorr_Shunt = full_autocorr_Shunt / full_autocorr_Shunt[mid_point_Shunt]

mid_point_Uniform = len(full_autocorr_Uniform) // 2
normalized_autocorr_Uniform = full_autocorr_Uniform / full_autocorr_Uniform[mid_point_Uniform]

# Compute lags for plotting
lags = np.arange(-len(hist_counts0) + 1, len(hist_counts0))

# Plot the normalized auto-correlation
plt.figure()
plt.plot(lags, normalized_autocorr_Hyper,'b', label = 'II = 0.3, GJ = True')
plt.plot(lags, normalized_autocorr_Shunt,'r',label = 'II = 0.0, GJ = False')
plt.plot(lags, normalized_autocorr_Uniform,'g',markersize = 4, label = 'II = 0.0, GJ = True')

plt.legend()
plt.title('Normalized Auto-correlation of Spike Histogram')
plt.xlabel('Lag')
plt.ylabel('Normalized Auto-correlation')
plt.ylim(-0.5, 1)
plt.savefig(fileID+'ING_Auto_Corr.png', format='png')
plt.savefig(fileID+'ING_Auto_Corr.eps', format='eps')
plt.show()


