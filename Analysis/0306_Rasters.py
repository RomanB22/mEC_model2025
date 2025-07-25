import numpy as np
from funcsAux import *
import matplotlib.pyplot as plt
import glob, os
from netpyne.analysis.tools import loadData
from itertools import compress
import matplotlib.gridspec as gridspec

spktFSAux = {}
i=0
for file in glob.glob("../0306_Pastoll/Ant&Ca_0_Kv3_1.0_Kv7_1.0_gsinInh7.0_gsinExc0.0_*data.pkl"):
    fileInfo = loadData(file)
    tSim = fileInfo['simData']['t']
    print(fileInfo['simConfig']['filename'][8:])
    spkt = np.array(fileInfo['simData']['spkt'])#[0:-1:10])
    spkid = np.array(fileInfo['simData']['spkid'])#[0:-1:10])
    fileID = '0306_Plots/ING_'

    # Time window of interest (first 250 ms)
    t_start = 0
    t_end = 250

    # Create a mask for the time window
    time_mask = (spkt >= t_start) & (spkt <= t_end)

    # Apply the time mask to spike times and IDs
    spkt_window = spkt[time_mask]
    spkid_window = spkid[time_mask]


    maskFS = spkid_window < 100 #mask for FS cells
    maskSC = spkid_window >= 100 #mask for SC cells

    # Subset spike times and neuron IDs based on the mask
    spkt_FS = spkt_window[maskFS]
    spkid_FS = spkid_window[maskFS]

    spkt_SC = spkt_window[maskSC]
    spkid_SC = spkid_window[maskSC]

    # --- FS Cell Rasters ---
    n_neurons_per_plot = 20
    for i in range(0, 100, n_neurons_per_plot):  # Iterate through FS cells (0-99)
        start_id = i
        end_id = min(i + n_neurons_per_plot, 100)  # Cap at 100
        neuron_mask = (spkid_FS >= start_id) & (spkid_FS < end_id)
        spkt_plot = spkt_FS[neuron_mask]
        spkid_plot = spkid_FS[neuron_mask]

        sorted_indices = np.argsort(spkid_plot)
        spkt_sorted = spkt_plot[sorted_indices]
        spkid_sorted = spkid_plot[sorted_indices]

        if len(spkt_sorted) > 0: #check if there are any spikes to plot. If not, the figure will be empty
            plt.figure()
            plt.plot(spkt_sorted, spkid_sorted, '.b', markersize=4)
            plt.xlabel('Time (ms)')
            plt.ylabel('Neuron ID (FS)')
            plt.title(f'Raster Plot of FS Neurons ({start_id}-{end_id-1})')
            plt.gca().set_yticks(np.arange(start_id, end_id, 1)) #set y ticks to be integers and in steps of 5
            plt.gca().set_yticklabels(np.arange(start_id, end_id, 1)) #set y tick labels to be integers and in steps of 5

            plt.tight_layout()  # Further adjust layout for tightness
            plt.savefig(fileID+str(start_id)+'to'+str(end_id)+'_RasterFS.png', format='png')
            plt.savefig(fileID+str(start_id)+'to'+str(end_id)+'_RasterFS.eps', format='eps')
            #plt.show()

    # --- SC Cell Rasters ---
    for i in range(100, 500, n_neurons_per_plot):  # Iterate through SC cells (100-499)
        start_id = i
        end_id = min(i + n_neurons_per_plot, 500)  # Cap at 500
        neuron_mask = (spkid_SC >= start_id) & (spkid_SC < end_id)
        spkt_plot = spkt_SC[neuron_mask]
        spkid_plot = spkid_SC[neuron_mask]

        sorted_indices = np.argsort(spkid_plot)
        spkt_sorted = spkt_plot[sorted_indices]
        spkid_sorted = spkid_plot[sorted_indices]

        if len(spkt_sorted) > 0: #check if there are any spikes to plot. If not, the figure will be empty
            plt.figure()
            plt.plot(spkt_sorted, spkid_sorted, '.r', markersize=4)
            plt.xlabel('Time (ms)')
            plt.ylabel('Neuron ID (SC)')
            plt.title(f'Raster Plot of SC Neurons ({start_id}-{end_id-1})')
            plt.gca().set_yticks(np.arange(start_id, end_id, 1)) #set y ticks to be integers and in steps of 5
            plt.gca().set_yticklabels(np.arange(start_id, end_id, 1)) #set y tick labels to be integers and in steps of 5

            plt.tight_layout()  # Further adjust layout for tightness
            plt.savefig(fileID+str(start_id)+'to'+str(end_id)+'_RasterSC.png', format='png')
            plt.savefig(fileID+str(start_id)+'to'+str(end_id)+'_RasterSC.eps', format='eps')
            #plt.show()
