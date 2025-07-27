"""
init.py

Starting script to run NetPyNE-based mEC model.
"""
import numpy as np
from netpyne import sim

###############################################################################
# cfg, netParams = sim.loadFromIndexFile('index.npjson') 
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='spatialModel/cfg.py', netParamsDefault='spatialModel/netParams.py')
sim.create(netParams = netParams, simConfig = cfg)
sim.simulate()
sim.analyze()
