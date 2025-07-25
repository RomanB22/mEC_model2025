"""
defs.py

Auxiliar file with some cell definitions and functions used in the mEC model.
"""
from netpyne import specs
import numpy as np
import gc
import os

def PVCell(cfg, gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams):
    netParamsAux = specs.NetParams()
    #To get a length from net capacitance most literature assumes neurons have 1uf/cm^2. This capacitance is in nF
    SAOrig=CapsOrig*1e-3
    SAGapJunct=CapsMod*1e-3
    #neuron does length in um so lets make it um^2
    SAumOrig=SAOrig*1e8
    SAumGapJunct=SAGapJunct*1e8
    #we're going to assume a diameter of 20um to get a radius of 10um
    #SA of a cyclinder(not including ends)=2*pi*r*L
    lengthsOrig=SAumOrig/(20*np.pi)
    lengthsGapJunct=SAumGapJunct/(20*np.pi)
    # Active conductances
    # Heterogeneous peak conductances for active currents. This gave me nanosiemens. I need siemens/cm^2 for it's conductance, so I need to multiple by e-9 to get nanosiemens then divide by the surface area in cm^2
    gLsSA = gLs*1e-9/SAOrig
    ConductWithGapJunctSA = ConductWithGapJunct*1e-9/SAGapJunct
    gNasSA = gNas*1e-9/SAOrig
    gNasWithGapJunctSA = gNas*1e-9/SAGapJunct
    gKv3sSA = gKv3s*1e-9/SAOrig
    gKv3sWithGapJunctSA = gKv3s*1e-9/SAGapJunct
    gKv7sSA = gKv7s*1e-9/SAOrig
    gKv7sWithGapJunctSA = gKv7s*1e-9/SAGapJunct
    # Fast Spiking PV+ Basket Cell from mEC
    for k in range(cfg.NumModels):
        cellRule = {'conds': {'cellType': 'FS'}, 'diversityFraction': 1./cfg.NumModels , 'secs': {}}  # cell rule dict        
        if cfg.OPTODRIVE==False:
            cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}} 
        else:
            cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}, 'pointps': {}}
            drive = cfg.g_sin if cfg.HETERDRVE==False else cfg.drives[k]
            cellRule['secs']['soma']['pointps']['optodrive'] = {'mod' : 'optodrive', 'gsin': drive, 'Ese': 0 , 'f': cfg.fsin}
        if cfg.NOISE==True:
            cellRule['secs']['soma']['pointps']['InVivoNoise'] = {'mod' : 'Gfluct','g_e0': cfg.MeanENoise, 'g_i0': cfg.MeanINoise, 'std_e': cfg.StdENoise, 'std_i': cfg.StdINoise}
        if cfg.GAP==True:
            cellRule['secs']['soma']['geom'] = {'diam': 20, 'L': lengthsGapJunct[k], 'cm': 1}   # soma geometry
            cellRule['secs']['soma']['mechs']['pas'] = {'g': ConductWithGapJunctSA[k], 'e': ReversPotWithGapJunct[k]}
            cellRule['secs']['soma']['mechs']['naG'] = {'gbar': gNasWithGapJunctSA[k], 'thm1': thm1s[k], 'thh2': thh2s[k] }
            cellRule['secs']['soma']['mechs']['kv7'] = {'gbar': gKv7sWithGapJunctSA[k], 'tha1': tha1s[k], 'ka1': SharedParams[13], 'ka2':SharedParams[14] }
            cellRule['secs']['soma']['mechs']['kv3'] = {'gbar': gKv3sWithGapJunctSA[k], 'thn1': thn1s[k], 'kn1':SharedParams[10], 'kn2':SharedParams[11] }
        else:
            cellRule['secs']['soma']['geom'] = {'diam': 20, 'L': lengthsOrig[k], 'cm': 1}   # soma geometry
            cellRule['secs']['soma']['mechs']['pas'] = {'g': gLsSA[k], 'e': ELs[k]}
            cellRule['secs']['soma']['mechs']['naG'] = {'gbar': gNasSA[k], 'thm1': thm1s[k], 'thh2': thh2s[k] }
            cellRule['secs']['soma']['mechs']['kv7'] = {'gbar': gKv7sSA[k], 'tha1': tha1s[k], 'ka1': SharedParams[13],'ka2':SharedParams[14]}
            cellRule['secs']['soma']['mechs']['kv3'] = {'gbar': gKv3sSA[k], 'thn1': thn1s[k], 'kn1':SharedParams[10], 'kn2':SharedParams[11] }    
        cellRule['secs']['soma']['vinit'] = np.random.uniform(-75,-65) # set initial membrane potential
        netParamsAux.cellParams['FS_'+str(k)] = cellRule
    gc.collect()
    return netParamsAux.cellParams

def SCell_HH(cfg):
    netParamsAux = specs.NetParams()
    # Stellate Cell (Modified Hodgkin Huxley)
    SCcell = {'secs': {}}
    if cfg.OPTODRIVE==False:
        SCcell['secs']['soma'] = {'geom': {}, 'mechs': {}} 
    else:
        SCcell['secs']['soma'] = {'geom': {}, 'mechs': {}, 'pointps': {}}
        SCcell['secs']['soma']['pointps']['optodrive'] = {'mod' : 'optodrive', 'gsin': cfg.g_sinExc, 'Ese': 0, 'f': cfg.fsin }
    if cfg.NOISE==True:
        SCcell['secs']['soma']['pointps']['InVivoNoise'] = {'mod' : 'Gfluct', 'g_e0': cfg.MeanENoise, 'g_i0': cfg.MeanINoise, 'std_e': cfg.StdENoise, 'std_i': cfg.StdINoise}
    SCcell['secs']['soma']['geom'] = {'diam': 18.2, 'L': 18.2, 'Ra': 150, 'cm':1}                           # soma geometry
    SCcell['secs']['soma']['mechs']['hh'] = {'gnabar': '0.12*normal(1,3e-2)', 'gkbar': '0.036*normal(1,3e-2)', 'gl': '0.0000357*normal(1.2,1e-1)', 'el': -68}  
    SCcell['secs']['soma']['vinit'] = np.random.uniform(-65,-58) # set initial membrane potential
    netParamsAux.cellParams['SC'] = SCcell

    gc.collect()
    return netParamsAux.cellParams

def SC_Mittal(cwd, cfg):
    netParamsAux = specs.NetParams()

    cellRule = netParamsAux.importCellParams(label='SC', conds={'cellType': 'SC'},
    fileName=cwd+'/cells/Stochastic_Bifurcations_Mittal_Narayanan/Figure3-6_SCmodel/Additive_Noise/Noise_Osc_july.hoc',
    cellName='SimpleNeuron', cellInstance = True)

    if cfg.OPTODRIVE==True:
        cellRule['secs']['soma']['pointps'] = {}
        cellRule['secs']['soma']['pointps']['optodrive'] = {'mod' : 'optodrive', 'gsin': cfg.g_sinExc, 'Ese': 0, 'f': cfg.fsin }

    del cellRule['secLists']

    gc.collect()
    return netParamsAux.cellParams
