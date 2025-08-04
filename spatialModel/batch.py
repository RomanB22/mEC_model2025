"""
batch.py 

Batch simulation for M1 model using NetPyNE

Contributors: romanbaravalle@gmail.com
"""
from netpyne.batch import Batch
from netpyne import specs
import os
import numpy as np
import random
random.seed(7894)

def runNetworks(NumNetworks=2):
    params = specs.ODict()
    # PV connectivity seeds
    seeds = [random.randrange(7894) for i in range(NumNetworks-1)]
    seeds.insert(0, 7894)  # First seed is always the same
    params[('seeds', 'brian2')] = seeds
    b = Batch(params=params, netParamsFile='spatialModel/netParams.py', cfgFile='spatialModel/cfg.py')
    b.method = 'grid'

    return b

def runDifferentReversalPot():
    params = specs.ODict()
    params['SYNAPSES'] = ['Hyper','Shunt','Uniform'] 
    b = Batch(params=params, netParamsFile='spatialModel/netParams.py', cfgFile='spatialModel/cfg.py')
    b.method = 'grid'

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin'):
    if type=='mpi_bulletin' or type=='mpi':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'spatialModel/init.py', 
            'skip': True}
    
# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    # b = runNetworks(NumNetworks=1)
    b = runDifferentReversalPot()
    b.batchLabel = 'batch_mEC_0'  # Label for the batch
    os.makedirs('./outputSpatial', exist_ok=True)
    b.saveFolder = './outputSpatial/'+b.batchLabel
    setRunCfg(b, 'mpi_bulletin')
    b.run() # run batch