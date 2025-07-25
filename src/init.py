"""
init.py

Starting script to run NetPyNE-based mEC model.
"""
import numpy as np
from netpyne import sim

###############################################################################
# cfg, netParams = sim.loadFromIndexFile('index.npjson') 
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='src/cfg.py', netParamsDefault='src/netParams.py')
sim.create(netParams = netParams, simConfig = cfg)
sim.simulate()
sim.analyze()
