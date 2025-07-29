import Inet.CreateNetworkParameters as Inet
import src.defs as defs
import numpy as np
from netpyne import specs
import os
try:
    from __main__ import cfg
except:
    from spatialModel.cfg import cfg

cwd = os.getcwd()  # Get current working directory

# Default network parameters
netParams = specs.NetParams()  # object of class NetParams to store the network parameters
netParams.defaultDelay = 0 #This is because it creates gap junctions with defaultDelay = 1 ms if not 
netParams.defaultThreshold = -30.0

netParams.shape = 'cuboid'
netParams.sizeX = cfg.xlength # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.ylength # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.zlength # z-dimension (horizontal depth) size in um


###############################################################################
## Create and load parameters for the PV network
###############################################################################
# gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, syns, delays, gms, synsgj, ggs = Inet.NetworkParams(NumNeurons=cfg.NPV, FactorTau=cfg.FactorTau, FactorKv3=cfg.FactorKv3, FactorKv7=cfg.FactorKv7, GapJunctProb=cfg.GapJunctProb, ChemycalConnProb=cfg.ChemycalConnProb, delaymin=.6, delaymax=1., meangms=0., sigmagms=1.,homogeneous=cfg.HOMOGENEOUS_PV, randomseed = cfg.seeds['brian2'])

gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, synsgj, ggs = Inet.ViaModelParams(NumNeurons=cfg.NPV, FactorTau=cfg.FactorTau, FactorKv3=cfg.FactorKv3, FactorKv7=cfg.FactorKv7, GapJunctProb=0,homogeneous=cfg.HOMOGENEOUS_PV, randomseed = cfg.seeds['brian2'])
# syns, delays, gms = Inet.NetConnectivity(ChemycalConnProb=cfg.ChemycalConnProb, delaymin=.6, delaymax=1., meangms=0., sigmagms=1.)
netParams.ConductOrig = gLs  # Original conductance for the PV+ cells
netParams.ReversPotOrig = ELs  # Original reversal potential for the PV+ cells
###############################################################################
## Cell types
###############################################################################
cellParamsPV = defs.PVCell(cfg, gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams)
cellParamsSC = defs.SCell_HH(cfg)
cellParamsSC_Mittal = defs.SC_Mittal(cwd, cfg)

netParams.cellParams = cellParamsPV | cellParamsSC_Mittal if cfg.Mittal else cellParamsPV | cellParamsSC # Combine all dictionaries

###############################################################################
# NETWORK PARAMETERS
###############################################################################
# Population parameters
netParams.popParams['FS'] = {'cellType': 'FS', 'diversity': True, 'density': cfg.PVdensity} # add dict with params for this pop
netParams.popParams['SC'] = {'cellType': 'SC', 'diversity': True, 'density': cfg.SCdensity} # add dict with params for this pop
#netParams.popParams['PYR'] = {'cellType': 'PYR', 'numCells': 1} 

###############################################################################
## VoltageClamp
###############################################################################
if cfg.Clamp==True:
	netParams.stimSourceParams['Vclamp'] = {'type': 'SEClamp', 'dur1': 1e9, 'amp1': cfg.Vclamp, 'rs': 1e-4}
    # Stimulation mapping parameters
	netParams.stimTargetParams['Vclamp->Cells'] = {
        'source': 'Vclamp',
        'sec': 'soma',
        'loc': 0.5,
        'conds': {'cellList': cfg.ClampCells}}

#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
    for key in [k for k in dir(cfg) if k.startswith('IClamp')]:
        params = getattr(cfg, key, None)
        [pop,sec,loc,start,dur,amp] = [params[s] for s in ['pop','sec','loc','start','dur','amp']]

        #cfg.analysis['plotTraces']['include'].append((pop,0))  # record that pop

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'IClamp', 'delay': start, 'dur': dur, 'amp': amp}
        
        # connect stim source to target
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key, 
            'conds': {'pop': pop},
            'sec': sec, 
            'loc': loc}

###############################################################################
## Synaptic mechs
###############################################################################
## Inhibitory synapses FS-> FS
tau_rise=0.3
tau_fall=2.
c_fall = 1/tau_fall
c_rise = 1/tau_rise
f = 1/( np.exp(-c_fall*np.log(c_rise/c_fall)/(c_rise-c_fall)) - np.exp(-c_rise*np.log(c_rise/c_fall)/(c_rise-c_fall)) )
U_SE=0.3; tau_d=100.
# Let's add synapses  
if cfg.WODepression==True:
    netParams.synMechParams['inhFSFS'] = {'mod': 'synact', 'tau_rise': tau_rise, 'tau_fall': tau_fall, 'f' : f, 'Es': cfg.Esyn_inh[cfg.SYNAPSES]}        
else:
    # gms = [i/U_SE for i in gms]
    netParams.synMechParams['inhFSFS'] = {'mod': 'synactdep', 'tau_rise': tau_rise, 'tau_fall': tau_fall, 'f' : f, 'xs' : 1, 'U_SE' : U_SE, 'tau_d' : tau_d, 'Es': cfg.Esyn_inh[cfg.SYNAPSES]}
# gms = [cfg.gmsScale*i for i in gms]
# ggs = [cfg.ggsScale*i for i in ggs]
# Connectivity parameters
netParams.connParams['FS->FS_chem'] = {
        'preConds': {'pop': 'FS'},         # presynaptic conditions
        'postConds': {'pop': 'FS'},        # postsynaptic conditions
        'sec':'soma',
        'probability': cfg.I2IProbchem,
        'weight': cfg.weightI2Ichem,                 
        'synMech': 'inhFSFS',                   # target inh synapse
         'delay': cfg.delaysI2Ichem}                    # delay

###############################################################################
## Inhibitory synapses FS-> SC

tau_riseExc=0.4
tau_fallExc=6.
c_fall = 1./tau_fallExc; c_rise = 1./tau_riseExc
norm_synExc = 1./( np.exp(-c_fall*np.log(c_rise/c_fall)/(c_rise-c_fall)) - np.exp(-c_rise*np.log(c_rise/c_fall)/(c_rise-c_fall)) )
netParams.synMechParams['inhFSSC'] = {'mod': 'synactdep', 'tau_rise': tau_riseExc, 'tau_fall': tau_fallExc, 'f' : norm_synExc, 'xs' : 1, 'U_SE' : U_SE, 'tau_d' : tau_d, 'Es': -65.}
# Connectivity parameters
netParams.connParams['FS->SC'] = {
        'preConds': {'pop': 'FS'},         # presynaptic conditions
        'postConds': {'pop': 'SC'},        # postsynaptic conditions
        'sec':'soma',
        'probability': cfg.I2EProb,
        #I'm subbing weight for conductance, gms was in nanosiemens and needs to be converted to uS
        'weight': cfg.weightI2E,                      # weight of each connection. Take into account that numpy and NEURON arguments are different (numpy args are mean and std for subjacent normal distribution, not the lognorm as in NEURON)
        'synMech': 'inhFSSC',                   # target inh synapse
        'delay': cfg.delaysI2E}                    # delay

###############################################################################
## Excitatory synapses SC -> FS

netParams.synMechParams['AMPA'] = {'mod': 'ExpSyn', 'tau': 1, 'e': 0}  # excitatory synaptic mechanism
# Synaptic mechanism parameters
netParams.synMechParams['NMDA'] = {'mod': 'Exp2Syn', 'tau1': 0.1, 'tau2': 5.0, 'e': 0}  # NMDA synaptic mechanism (NOT IMPLEMENTED YET)
# Connectivity parameters
netParams.connParams['SC->FS'] = {
        'preConds': {'pop': 'SC'},         # presynaptic conditions
        'postConds': {'pop': 'FS'},        # postsynaptic conditions
        'sec':'soma',
        'probability': cfg.I2EProb,
        #am subbing weight for conductance, gms was in nanosiemens and needs to be converted to uS
        'weight': cfg.weightI2E,              # weight of each connection
        'synMech': 'AMPA',                   # target exc synapse
         'delay': cfg.delaysI2E}                   # delay