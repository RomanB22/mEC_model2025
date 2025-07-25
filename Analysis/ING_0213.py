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
maxPower= 900 # 100. # Sets the maximum for the scalogram. Varies according to the power of the signal. Change it for better visualization
colorbar=True
numBins = 2*125+1
###################################
fs = np.arange(minFreq,maxFreq,3.)
levelsWav = np.linspace(0., maxPower, 20, endpoint=True)


i = 0
spktFSAux = {}
spktSCAux = {}
FScellTrace={}
SCcellTrace={}
for file in glob.glob("../0213_Higher_E_Drive/Ant&Ca_0_Kv3_1.0_Kv7_1.0_gsinInh7.0_gsinExc0.0_*Hyper_Noise_F_data.pkl"):
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][8:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    FScellTrace['0'] = fileInfo['simData']['V_soma']['cell_70']
    FScellTrace['1'] = fileInfo['simData']['V_soma']['cell_80']
    FScellTrace['2'] = fileInfo['simData']['V_soma']['cell_45']
    FScellTrace['3'] = fileInfo['simData']['V_soma']['cell_95']
    FScellTrace['4'] = fileInfo['simData']['V_soma']['cell_15']
    FScellTrace['5'] = fileInfo['simData']['V_soma']['cell_35']
    maskFS = np.array(fileInfo['simData']['spkid'])< 100
    spktFSAux['0'] = np.array( list( compress(fileInfo['simData']['spkt'], maskFS) ) )
    tSim = fileInfo['simData']['t']
    i += 1

spktFS0 = spktFSAux['0']
cellTrFS0 = FScellTrace['0']
cellTrFS1 = FScellTrace['1']
cellTrFS2 = FScellTrace['2']
cellTrFS3 = FScellTrace['3']
cellTrFS4 = FScellTrace['4']
cellTrFS5 = FScellTrace['5']

fileID = '0213_Plots/0213_ING_7nS_'

# plot it
fig = plt.figure(figsize=(12, 6))
outer = gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.1)

inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.1, height_ratios=[1,3])
 
ax0 = plt.Subplot(fig, inner[0])
t = np.linspace(sim_time-2.*125.,sim_time,2*125+1)
ax0.plot(t, gsin*np.sin(2*np.pi*f*t*1e-3-np.pi/2),color = 'k')
ax0.set_xlim(sim_time-2.*125.,sim_time)
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.set_xticks([])
ax0.set_yticks([])
fig.add_subplot(ax0)
   
ax1 = plt.Subplot(fig, inner[1])
histFS0, foo1, foo = ax1.hist(spktFS0, np.linspace(sim_time-2.*125.,sim_time, numBins),color = 'k', label='FS Stim Only')
ax1.set_xlim(sim_time-2.*125.,sim_time); tmp = ax1.set_ylim(0.,20.)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.spines['left'].set_visible(True)
ax1.set_xticks([])
ax1.set_xticklabels([])
ax1.set_ylabel('Counts')
ax1.legend(loc='upper right')
fig.add_subplot(ax1)

# WAVELET PLOTS
inner1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[1], wspace=0.3, hspace=0.1)

ax4 = plt.Subplot(fig, inner1[0])
cwtPowFS, cwtPhaseFS = compPWT(histFS0,1e-3,fs=fs)
plotPow(cwtPowFS,dt=1e-3,fs=fs,levels=levelsWav, cmap="hot", ax=ax4, fig=fig, colorbar=colorbar)
#plotPow(cwtPowFS,dt=dt_rec,fs=fs, cmap="hot", ax=ax3, fig=fig, colorbar=colorbar)
ax4.set_xticks([]) 
ax4.set_ylabel("Freq (Hz)")
#ax4.set_ylim(50,180)
fig.add_subplot(ax4)


fig.tight_layout()
fig.savefig(fileID+'Hist.png', bbox_inches='tight')
fig.savefig(fileID+'Hist.eps', bbox_inches='tight')

