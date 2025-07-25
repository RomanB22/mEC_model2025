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
maxFreq=350.
maxPower= 500 # 100. # Sets the maximum for the scalogram. Varies according to the power of the signal. Change it for better visualization
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
for file in glob.glob("../0226_Pastoll/Ant*gsinInh7.0*gsinExc0.0*Hyper*_Noise_F*.pkl"):
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][8:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    FScellTrace['0'] = fileInfo['simData']['V_soma']['cell_10']
    FScellTrace['1'] = fileInfo['simData']['V_soma']['cell_45']
    FScellTrace['2'] = fileInfo['simData']['V_soma']['cell_55']
    FScellTrace['3'] = fileInfo['simData']['V_soma']['cell_95']
    maskFS = np.array(fileInfo['simData']['spkid'])< 100
    spktFSAux['0'] = np.array( list( compress(fileInfo['simData']['spkt'], maskFS) ) )
    i += 1

fileID = '0227_Plots/0227_Pastoll_Const_11_'

spktFS0 = spktFSAux['0']
cellTrFS0 = FScellTrace['0']
cellTrFS1 = FScellTrace['1']
cellTrFS2 = FScellTrace['2']
cellTrFS3 = FScellTrace['3']

#print(type(cellTrFS0), len(cellTrFS0), type(FScellTrace['0']))


for file in glob.glob("../0226_Pastoll/Const_11_Ant&Ca_0_Kv3_1.0_Kv7_1.0_gsinInh7.0_gsinExc3.0_EI_0.3_IE_0.2_II_0.0_GAP_T_*_data.pkl"):
    fileInfo = loadData(file)
    fileName = fileInfo['simConfig']['filename'][8:]
    print(fileName)
    sim_time = fileInfo['simConfig']['duration']
    FScellTrace['4'] = fileInfo['simData']['V_soma']['cell_10']
    FScellTrace['5'] = fileInfo['simData']['V_soma']['cell_45']
    FScellTrace['6'] = fileInfo['simData']['V_soma']['cell_55']
    FScellTrace['7'] = fileInfo['simData']['V_soma']['cell_95']
    maskFS = np.array(fileInfo['simData']['spkid'])< 100
    spktFSAux['1'] = np.array( list( compress(fileInfo['simData']['spkt'], maskFS) ) )
    SCcellTrace['0'] = fileInfo['simData']['V_soma']['cell_120']# [1,20,40,50,90,105,110,130,230,250,270,290]
    SCcellTrace['1'] = fileInfo['simData']['V_soma']['cell_140']
    SCcellTrace['2'] = fileInfo['simData']['V_soma']['cell_225']
    SCcellTrace['3'] = fileInfo['simData']['V_soma']['cell_245']
    SCcellTrace['4'] = fileInfo['simData']['V_soma']['cell_275']
    #SCcellTrace['5'] = fileInfo['simData']['V_soma']['cell_400']
    maskSC = np.array(fileInfo['simData']['spkid'])>= 100
    tSim = fileInfo['simData']['t']
    spktSCAux['0'] = np.array( list( compress(fileInfo['simData']['spkt'], maskSC) ) )
    i += 1


spktFS1 = spktFSAux['1']
cellTrFS4 = FScellTrace['4']
cellTrFS5 = FScellTrace['5']
cellTrFS6 = FScellTrace['6']
cellTrFS7 = FScellTrace['7']

#sim.allSimData.LFP[:,0]
spktSC = spktSCAux['0']
cellTrSC0 = SCcellTrace['0']
cellTrSC1 = SCcellTrace['1']
cellTrSC2 = SCcellTrace['2']
cellTrSC3 = SCcellTrace['3']

# plot it
fig = plt.figure(figsize=(10, 8))
outer = gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.1)

inner = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.1, height_ratios=[1,4,4,4])
 
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
histFS0, foo1, foo = ax1.hist(spktFS0, np.linspace(sim_time-2.*125.,sim_time, numBins),color = 'k', label='FS. Stimulate FS only')
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

ax2 = plt.Subplot(fig, inner[2])
histFS1, foo, foo = ax2.hist(spktFS1, np.linspace(sim_time-2.*125.,sim_time, numBins), color = 'k', label='FS. Stimulate FS and SC')
ax2.set_xlim(sim_time-2.*125.,sim_time); tmp = ax2.set_ylim(0.,20.)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(True)
ax2.set_ylabel('Counts')
ax2.set_xticks([])
ax2.set_xticklabels([])
ax2.legend(loc='upper right')
fig.add_subplot(ax2)

ax3 = plt.Subplot(fig, inner[3])
histSC, foo, foo = ax3.hist(spktSC, np.linspace(sim_time-2.*125.,sim_time, numBins), color = 'k', label='SC. Stimulate FS and SC')
ax3.set_xlim(sim_time-2.*125.,sim_time); tmp = ax3.set_ylim(0.,30.)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(True)
ax3.set_ylabel('Counts')
ax3.set_xlabel('t (ms)')
ax3.legend(loc='upper right')

fig.add_subplot(ax3)

# WAVELET PLOTS
inner1 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[1], wspace=0.3, hspace=0.1)

ax4 = plt.Subplot(fig, inner1[0])
cwtPowFS0, cwtPhaseFS0 = compPWT(histFS0,1e-3,fs=fs)

maxPowerFS0 = np.max(np.max(cwtPowFS0))
print('Maximum Power in the FS0 case',maxPowerFS0)

# Find the indices of the maximum power
max_power_indexFS0 = np.unravel_index(np.argmax(cwtPowFS0), cwtPowFS0.shape)
max_freq_indexFS0, max_time_indexFS0 = max_power_indexFS0

# Get the frequency and time at the maximum power
max_freqFS0 = fs[max_freq_indexFS0]
print(f"Maximum power for FS firing occurs for ING at frequency: {max_freqFS0} Hz")

plotPow(cwtPowFS0,dt=1e-3,fs=fs,levels=levelsWav, cmap="hot", ax=ax4, fig=fig, colorbar=colorbar)
#plotPow(cwtPowFS,dt=dt_rec,fs=fs, cmap="hot", ax=ax3, fig=fig, colorbar=colorbar)
ax4.set_xticks([]) 
ax4.set_ylabel("Freq (Hz)")
fig.add_subplot(ax4)

ax5 = plt.Subplot(fig, inner1[1])
cwtPowFS1, cwtPhaseFS1 = compPWT(histFS1,1e-3,fs=fs)

maxPowerFS1 = np.max(np.max(cwtPowFS1))
print('Maximum Power in the FS1 case',maxPowerFS1)

# Find the indices of the maximum power
max_power_indexFS1 = np.unravel_index(np.argmax(cwtPowFS1), cwtPowFS1.shape)
max_freq_indexFS1, max_time_indexFS1 = max_power_indexFS1

# Get the frequency and time at the maximum power
max_freqFS1 = fs[max_freq_indexFS1]
print(f"Maximum power for FS firing occurs for PING at frequency: {max_freqFS1} Hz")
plotPow(cwtPowFS1,dt=1e-3,fs=fs,levels=levelsWav, cmap="hot", ax=ax5, fig=fig, colorbar=colorbar)
#plotPow(cwtPowSC,dt=dt_rec,fs=fs, cmap="hot", ax=ax4, fig=fig, colorbar=colorbar)
ax5.set_xticks([]) 
ax5.set_ylabel("Freq (Hz)")    
fig.add_subplot(ax5)

ax6 = plt.Subplot(fig, inner1[2])
cwtPowSC, cwtPhaseSC = compPWT(histSC,1e-3,fs=fs)

maxPowerSC = np.max(np.max(cwtPowSC))
print('Maximum Power in the SC case',maxPowerSC)

# Find the indices of the maximum power
max_power_indexSC = np.unravel_index(np.argmax(cwtPowSC), cwtPowSC.shape)
max_freq_indexSC, max_time_indexSC = max_power_indexSC

# Get the frequency and time at the maximum power
max_freqSC = fs[max_freq_indexSC]
print(f"Maximum power for SC occurs for PING at frequency: {max_freqSC} Hz")
plotPow(cwtPowSC,dt=1e-3,fs=fs,levels=levelsWav, cmap="hot", ax=ax6, fig=fig, colorbar=colorbar)
#plotPow(cwtPowSC,dt=dt_rec,fs=fs, cmap="hot", ax=ax4, fig=fig, colorbar=colorbar)
ax6.set_xticks([]) 
ax6.set_ylabel("Freq (Hz)")    
fig.add_subplot(ax6)

fig.tight_layout()
fig.savefig(fileID+'Hist.png', bbox_inches='tight')
fig.savefig(fileID+'Hist.eps', bbox_inches='tight')
#fig.savefig(fileID+'.svg', bbox_inches='tight')


# PLOTTING TRACES 1
fig = plt.figure(figsize=(10, 8))
outer = gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.1)

inner = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.1, height_ratios=[1,4,4,4])


ax0 = plt.Subplot(fig, inner[0])
t = np.linspace(sim_time-4.*125.,sim_time,2*125+1)
ax0.plot(t, gsin*np.sin(2*np.pi*f*t*1e-3-np.pi/2),color = 'k')
ax0.set_xlim(sim_time-4.*125.,sim_time)
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.set_xticks([])
ax0.set_yticks([])
fig.add_subplot(ax0)

ax1 = plt.Subplot(fig, inner[1], label='FS - Stimulate only FS')
y = cellTrFS0
ax1.plot(tSim, y, 'b-', label='FS - Stimulate only FS')
ax1.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax1)

ax2 = plt.Subplot(fig, inner[2], label='FS - Stimulate FS and SC')
y = cellTrFS4
ax2.plot(tSim, y, 'b-', label='FS - Stimulate FS and SC')
ax2.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax2)

ax3 = plt.Subplot(fig, inner[3], label='SC - Stimulate FS and SC')
y = cellTrSC0
ax3.plot(tSim, y, 'b-', label='SC - Stimulate FS and SC')
ax3.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax3)

inner = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=outer[1], wspace=0.1, hspace=0.1, height_ratios=[1,4,4,4])

ax0 = plt.Subplot(fig, inner[0])
t = np.linspace(sim_time-4.*125.,sim_time,2*125+1)
ax0.plot(t, gsin*np.sin(2*np.pi*f*t*1e-3-np.pi/2),color = 'k')
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.set_xlim(sim_time-4.*125.,sim_time)
ax0.set_xticks([])
ax0.set_yticks([])
fig.add_subplot(ax0)

ax1 = plt.Subplot(fig, inner[1])
y = cellTrFS1
ax1.plot(tSim, y, 'b-', label='FS - Stimulate only FS')
ax1.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax1)


ax2 = plt.Subplot(fig, inner[2])
y = cellTrFS5
ax2.plot(tSim, y, 'b-', label='FS - Stimulate FS and SC')
ax2.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax2)


ax3 = plt.Subplot(fig, inner[3])
y = cellTrSC1
ax3.plot(tSim, y, 'b-', label='SC - Stimulate FS and SC')
ax3.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax3)

fig.tight_layout()
#fileName = '0110_Plots/FS_gsin5nS_gsinEx2nS_Hyper_Traces'
fig.savefig(fileID+'Traces1.png', bbox_inches='tight')
fig.savefig(fileID+'Traces1.eps', bbox_inches='tight')
#fig.savefig(fileID+'.svg', bbox_inches='tight')

#Plotting Traces 2

fig = plt.figure(figsize=(10, 8))
outer = gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.1)

inner = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.1, height_ratios=[1,4,4,4])


ax0 = plt.Subplot(fig, inner[0])
t = np.linspace(sim_time-4.*125.,sim_time,2*125+1)
ax0.plot(t, gsin*np.sin(2*np.pi*f*t*1e-3-np.pi/2),color = 'k')
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.set_xticks([])
ax0.set_yticks([])
ax0.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax0)

ax1 = plt.Subplot(fig, inner[1], label='FS - Stimulate only FS')
y = cellTrFS2
ax1.plot(tSim, y, 'b-', label='FS - Stimulate only FS')
ax1.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax1)

ax2 = plt.Subplot(fig, inner[2], label='FS - Stimulate FS and SC')
y = cellTrFS6
ax2.plot(tSim, y, 'b-', label='FS - Stimulate FS and SC')
ax2.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax2)

ax3 = plt.Subplot(fig, inner[3], label='SC - Stimulate FS and SC')
y = cellTrSC2
ax3.plot(tSim, y, 'b-', label='SC - Stimulate FS and SC')
ax3.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax3)

inner = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=outer[1], wspace=0.1, hspace=0.1, height_ratios=[1,4,4,4])

ax0 = plt.Subplot(fig, inner[0])
t = np.linspace(sim_time-4.*125.,sim_time,2*125+1)
ax0.plot(t, gsin*np.sin(2*np.pi*f*t*1e-3-np.pi/2),color = 'k')
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.set_xlim(sim_time-4.*125.,sim_time)
ax0.set_xticks([])
ax0.set_yticks([])
fig.add_subplot(ax0)

ax1 = plt.Subplot(fig, inner[1])
y = cellTrFS3
ax1.plot(tSim, y, 'b-', label='FS - Stimulate only FS')
ax1.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax1)


ax2 = plt.Subplot(fig, inner[2])
y = cellTrFS7
ax2.plot(tSim, y, 'b-', label='FS - Stimulate FS and SC')
ax2.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax2)


ax3 = plt.Subplot(fig, inner[3])
y = cellTrSC3
ax3.plot(tSim, y, 'b-', label='SC - Stimulate FS and SC')
ax3.set_xlim(sim_time-4.*125.,sim_time)
fig.add_subplot(ax3)

fig.tight_layout()
#fileName = '0110_Plots/FS_gsin5nS_gsinEx2nS_Hyper_Traces_2'
fig.savefig(fileID+'Traces2.png', bbox_inches='tight')
fig.savefig(fileID+'Traces2.eps', bbox_inches='tight')
#fig.savefig(fileName+'.svg', bbox_inches='tight')

