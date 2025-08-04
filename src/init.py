"""
init.py

Starting script to run NetPyNE-based mEC model.
"""
import numpy as np
import src.defs as defs
from netpyne import sim

###############################################################################
# cfg, netParams = sim.loadFromIndexFile('index.npjson') 
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='src/cfg.py', netParamsDefault='src/netParams.py')
sim.create(netParams = netParams, simConfig = cfg)
sim.simulate()
# Calculate cell-to-cell connectivity matrix
connMatrix, pre, post = sim.analysis.network._plotConnCalculateFromSim(
    includePre=['all'],
    includePost=['all'],
    feature='numConns',
    orderBy='gid',
    groupBy='cell',
    groupByIntervalPre=None,
    groupByIntervalPost=None,
    synOrConn='conn',
    synMech=None,
    removeWeightNorm = False,
    logPlot=False,
)
sim.allSimData.ConnMatrix = connMatrix
sim.analyze()
if cfg.PlotWavelet:
    trange = [cfg.delayStim-100, cfg.delayStim+cfg.durationStim+100] # Change for theta drive
    defs.CalculateAndPlotWaveletClamp(cfg, sim, trange = trange, frange=[50, 300])

print('completed simulation...')