from netpyne import specs
import numpy as np

# Simulation options
cfg = specs.SimConfig()       # object of class SimConfig to store simulation configuration

###############################################################################
## Simulation parameters
###############################################################################
cfg.ThetaCycles = 4          # Number of theta cycles to simulate
cfg.Theta2Plot = 4          # Number of theta cycles to plot
cfg.duration = cfg.ThetaCycles*125.          # Duration of the simulation, in ms
cfg.dt = 1e-2                # Internal integration timestep to use
cfg.hParams = {'celsius': 23, 'v_init': -80}  
cfg.saveFolder = 'output'  # Folder to save output
cfg.simLabel = 'mEC_0'  # Simulation label, used in output file names
cfg.validateNetParams = False
cfg.verbose = False           # Show detailed messages
# cfg.progressBar = 0       # Show progress bar
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
cfg.saveCellConns = False
cfg.saveJson = True
cfg.savePickle = False

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------
cfg.addIClamp = False

cfg.IClamp1 = {'pop': 'FS', 'sec': 'soma', 'loc': 0.5, 'start': 50, 'dur': 100, 'amp': 0.50}
cfg.IClamp2 = {'pop': 'SC', 'sec': 'soma', 'loc': 0.5, 'start': 50, 'dur': 100, 'amp': 0.50}

###############################################################################
## SimParams
############################################################################### 
# PV+ cells properties 
cfg.NPV=100 # Number of Inhibitory neurons
cfg.HOMOGENEOUS_PV = False 
cfg.NumModelsPV = 1 if (cfg.HOMOGENEOUS_PV and not cfg.GAP) else cfg.NPV # With gap junctions and homogeneity I have an error in the code
cfg.FactorTau, cfg.FactorKv3, cfg.FactorKv7 = 1, 1, 1   # To modify the activation curves for ion channels and the membrane time constant

# Stellate cells properties 
cfg.Mittal = False # If True, uses the Mittal et al. model for Stellate cells
cfg.NSC=4*cfg.NPV # Number of Stellate cells
cfg.HOMOGENEOUS_SC = True 
cfg.NumModelsSC = 1 if cfg.HOMOGENEOUS_SC else 157 # Load all the valid SC models
cfg.SCidx = 0 # Which model to load if using homogeneous population

if cfg.Mittal==False: cfg.HOMOGENEOUS_SC=True

# Optogenetic drive                                                                                                                                                                                                                                                                                                                            
cfg.OPTODRIVE=True                                                   
cfg.g_sin = 7.*1e-3 # Optogenetic conductance for the inhibitory population                  
cfg.g_sinExc = 3.*1e-3 # Optogenetic conductance for the excitatory population.
cfg.fsin=8  # Optogenetic sinusoidal stimulation, in Hz
# Heterogeneous optogenetic drive for the PV+ cells
cfg.HETERDRIVE = True
np.random.seed(cfg.seeds['opto'])  # Fixed seed
cfg.drives = np.random.normal(loc=cfg.g_sin, scale=0.1*cfg.g_sin, size=cfg.NPV)

# Noisy stimulation
cfg.NOISE=False # If true, adds a random fluctuating current simulating in-vivo like noise. TODO: Adjust parameters for mEC region 
cfg.MeanENoise, cfg.MeanINoise, cfg.StdENoise, cfg.StdINoise = 0, 0*0.0001, 0.000001, 0.000001                                                                                                                                                                  

# Voltage clamp parameters
cfg.Clamp=False                                                                  
cfg.Vclamp = 0 # Voltage at which clamp, if cfg.Clamp==True                                                                                                       

# Synapses parameters                                                                                                                                                                                                                                 
cfg.gmsScale = 1 # Scaling for the synaptic conductances                                                                                                                                         
cfg.ggsScale = 1 # Scaling for gap junction conductances                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
cfg.GAP=True                                                                                                                                                                            
cfg.WODepression=False # Without short-term depression                                                                                                                                                                                                                                                                                                      
cfg.SYNAPSES = 'Hyper' # ['Hyper','Shunt','Uniform'] 
cfg.Esyn_inh = {'Hyper': -75., 'Shunt': -55., 'Uniform': 'uniform(-70,-55)'}
cfg.GapJunctProb, cfg.ChemycalConnProb = 0.05, 0.3 # Probability connections from FS PV+ to FS PV+                                                           
cfg.ConnProbIE, cfg.ConnProbEI = 0.2, 0.3 # FS PV+ to SC and SC to FS PV+ probability connections
cfg.Weight_E2I = 0.00003

###############################################################################
## Recording and analysis
###############################################################################
# The index of cells are: First 100 are PV, the rest are SC 
cfg.ClampCells = [1, 20, 130]#[1,20,40,50,90,105,110,120,130,230,250,270,290,300,350]
recordedCells2 = [1, 11, 120, 130] #[1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,350,400]
# recordedCells3 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,350,400]
for_raster = ['FS', 'SC'] # [5,11,12,15,18,25,29,35,43,50,55,56,57,65,67,68,70,71,72,76,78,81,93] # Cells to include in the raster plot
offset = 0
timeRange = [(cfg.ThetaCycles-cfg.Theta2Plot)*125.,cfg.duration+offset]
cfg.recordCells = recordedCells2
cfg.recordTraces = {'V_soma':{'sec':'soma','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.analysis['plotRaster'] = {'include': for_raster,'saveFig': True, 'timeRange': timeRange}                  # Plot a raster
cfg.analysis['plotSpikeHist'] = {'include': for_raster, 'saveFig': True, 'timeRange': timeRange, 'binSize': 1, 'measure': 'rate'}                  # Plot a Spike Histogram
# cfg.analysis['plotTraces'] = {'include': recordedCells2, 'saveFig': True, 'timeRange': timeRange}  # Plot recorded traces for this list of cells
