"""
init.py

Starting script to run NetPyNE-based mEC model.
"""
import numpy as np
from netpyne import sim
import src.defs as defs

###############################################################################
# cfg, netParams = sim.loadFromIndexFile('index.npjson') 
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='spatialModel/cfg.py', netParamsDefault='spatialModel/netParams.py')
sim.initialize(simConfig = cfg, netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations
if cfg.GAP:
    print('Creating gap junctions...')
    defs.GapJunctSpatialConnectivity(sim, netParams, gLmin=2., ggapHigh=1.3, full_dist=False, sigma=.4)
    print('Gap junctions created.')
sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.runSim()                      			# run parallel Neuron simulation  
sim.gatherData()                  			# gather spiking data and cell info from each node
sim.saveData()                    			# save params, cell info and sim output to file (pickle,mat,txt,etc)#
sim.analysis.plotData()         			# plot spike raster etc
