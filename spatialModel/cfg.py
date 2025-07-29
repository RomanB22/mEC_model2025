from netpyne import specs
import numpy as np

# Simulation options
cfg = specs.SimConfig()       # object of class SimConfig to store simulation configuration

###############################################################################
## Simulation parameters
###############################################################################
cfg.ThetaCycles = 12          # Number of theta cycles to simulate
cfg.Theta2Plot = 2          # Number of theta cycles to plot
cfg.duration = cfg.ThetaCycles*125.          # Duration of the simulation, in ms
cfg.dt = 1e-2                # Internal integration timestep to use
cfg.hParams = {'celsius': 23, 'v_init': -80}  
cfg.saveFolder = 'outputSpatial'  # Folder to save output
cfg.simLabel = 'mEC_0'  # Simulation label, used in output file names
cfg.validateNetParams = False
cfg.verbose = False           # Show detailed messages
# cfg.progressBar = 0       # Hide progress bar
cfg.recordStep = cfg.dt        # Step size in ms to save data (e.g. V traces, LFP, etc)
cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321, 'cell': 4321, 'brian2': 7894, 'opto': 42} # Random seeds for reproducibility. brian2 seed is for the PV network.
cfg.saveDataInclude = ['simData', 'simConfig', 'net']  # Which data to save in the output file
cfg.printRunTime = 0.1 # Print run time every 0.1 seconds
cfg.recordTime = False  
cfg.createNEURONObj = True
cfg.createPyStruct = True
cfg.backupCfgFile = None #['cfg.py', 'backupcfg/']
cfg.gatherOnlySimData = False
cfg.saveCellSecs = False
cfg.saveCellConns = True
cfg.saveJson = False
cfg.savePickle = True

###############################################################################
## SimParams
############################################################################### 
# Synapses parameters                                                                                                                                                                                                                                 
cfg.gmsScale = 1 # Scaling for the synaptic conductances                                                                                                                                         
cfg.ggsScale = 1 # Scaling for gap junction conductances                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
cfg.GAP=True # CHECK IT JUST IN CASE           
cfg.GapJunctMaxDist = 150 # Maximum distance for gap junction connections, in um                                                                                                                                                                
cfg.WODepression=False # Without short-term depression                                                                                                                                                                                                                                                                                                      
cfg.SYNAPSES = 'Hyper' # ['Hyper','Shunt','Uniform'] 
cfg.Esyn_inh = {'Hyper': -75., 'Shunt': -55., 'Uniform': 'uniform(-70,-55)'}
cfg.GapJunctProb, cfg.ChemycalConnProb = 0.01, 0.2 # Probability connections from FS PV+ to FS PV+  

cfg.Weight_E2I = 'lognormal(0.1*1.65, 0.1*2.17)*1e-3' # Weight from excitatory to inhibitory cells
cfg.weightI2E = 'lognormal(0.1*1.65, 0.1*2.17)*1e-3' # Weight from inhibitory to excitatory cells
cfg.weightI2Ichem = 'lognormal(0.1*1.65,2.17)*1e-3' # Weight from inhibitory to inhibitory cells

cfg.delaysE2I = '0.6+(1-0.6)*uniform(0,1)' # Delay from excitatory to inhibitory cells
cfg.delaysI2E = '0.6+(1-0.6)*uniform(0,1)' # Delay from inhibitory to excitatory cells
cfg.delaysI2Ichem = '0.6+(1-0.6)*uniform(0,1)'# Delay from inhibitory to inhibitory cells

# Could be 2D or 3D distance: dist_2D, dist_3D
cfg.E2IProb = 'exp(-(dist_3D**2)/(2*98.36**2))'
cfg.I2EProb = 'exp(-(dist_3D**2)/(2*98.36**2))'
cfg.I2IProbchem = 'exp(-(dist_3D**2)/(2*98.36**2))' 

# Geometry and density
cfg.xlength = 800 # x-dimension (horizontal length) size in um
cfg.ylength = 300 # y-dimension (vertical height or cortical depth) size in um
cfg.zlength = 800 # z-dimension (horizontal depth) size in um
cfg.PVdensity = 1672 # Density of PV+ cells in cells/mm^3. (Bjerke et al. 2020)
cfg.SCdensity = 4*cfg.PVdensity # Density of PV+ cells in cells/mm^3
volume = cfg.xlength * cfg.ylength * cfg.zlength / 1e9 # Volume in mm^3

# PV+ cells properties
cfg.NPV = int(cfg.PVdensity * volume)  # Number of Inhibitory neurons. NOT USED HERE
cfg.HOMOGENEOUS_PV = False 
cfg.NumModelsPV = 1 if (cfg.HOMOGENEOUS_PV and not cfg.GAP) else cfg.NPV # With gap junctions and homogeneity I have an error in the code
cfg.FactorTau, cfg.FactorKv3, cfg.FactorKv7 = 1, 1, 1   # To modify the activation curves for ion channels and the membrane time constant

# Stellate cells properties 
cfg.Mittal = True # If True, uses the Mittal et al. model for Stellate cells
cfg.NSC = int(cfg.SCdensity * volume)  # Number of Stellate cells. NOT USED HERE
cfg.HOMOGENEOUS_SC = False 
cfg.NumModelsSC = 1 if cfg.HOMOGENEOUS_SC else 157 # Load all the valid SC models
cfg.SCidx = 0 # Which model to load if using homogeneous population

# Optogenetic drive                                                                                                                                                                                                                                                                                                                            
cfg.OPTODRIVE=True                                                   
cfg.g_sin = 4.*1e-3 # Optogenetic conductance for the inhibitory population                  
cfg.g_sinExc = 6.*1e-3 # Optogenetic conductance for the excitatory population.
cfg.fsin=8  # Optogenetic sinusoidal stimulation, in Hz
# Heterogeneous optogenetic drive for the PV+ cells
cfg.HETERDRIVE = True
np.random.seed(cfg.seeds['opto'])  # Fixed seed
cfg.drives = np.random.normal(loc=cfg.g_sin, scale=0.1*cfg.g_sin, size=cfg.NPV)

# Noisy stimulation
cfg.NOISE=False # If true, adds a random fluctuating current simulating in-vivo like noise. TODO: Adjust parameters for mEC region 
cfg.MeanENoise, cfg.MeanINoise, cfg.StdENoise, cfg.StdINoise = 0, 0*0.0001, 0.000001, 0.000001                                                                                                                                                                  

# Voltage clamp parameters
cfg.Clamp=True                                                                  
cfg.Vclamp = 0 # Voltage at which clamp, if cfg.Clamp==True                                                                                                       

###############################################################################
## Recording and analysis
###############################################################################
# The index of cells are: First 100 are PV, the rest are SC 
cfg.ClampCells = [0, 50, 321, 400]
recordedCells2 = [1, 11, 120, 130, 500] 
indicesCell = [i for i in range(2)]
recordedCells3 = [(pop, i) for i in indicesCell for pop in ['FS', 'SC']] 
for_raster = ['FS', 'SC'] 
indicesFS = [i for i in range(20)]
indicesSC = [i for i in range(80)]
for2Ddistribution = [('FS', i) for i in indicesFS] # Cells to include in the 2D distribution plot
[for2Ddistribution.append(('SC', i)) for i in indicesSC]
offset = 0
timeRange = [(cfg.ThetaCycles-cfg.Theta2Plot)*125.,cfg.duration+offset]
cfg.recordCells = recordedCells2
cfg.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.analysis['plotRaster'] = {'include': for_raster,'saveFig': True, 'timeRange': timeRange}                  # Plot a raster
cfg.analysis['plotSpikeHist'] = {'include': ['FS', 'SC'], 'saveFig': True, 'timeRange': timeRange, 'binSize': 1, 'measure': 'rate'}                  # Plot a Spike Histogram
cfg.analysis['plotTraces'] = {'include': recordedCells3, 'saveFig': True, 'timeRange': timeRange}  # Plot recorded traces for this list of cells
cfg.analysis['plot2Dnet'] = {'include': for2Ddistribution, 'view': 'xz', 'saveFig': True, 'showFig': True, 'lineWidth': 0.1, 'showConns': True} # Plot 2D cells and connections