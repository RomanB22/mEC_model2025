import numpy as np
from scipy.stats import truncnorm
import random 
from brian2 import NeuronGroup, Synapses, seed, prefs, BrianLogger
import os

BrianLogger.log_level_error()

# See end of file for Usage

def gen_syns_samepop(N,p):
    net = NeuronGroup(N, "dv/dt = 0 : volt", threshold="v>1.", method="euler");
    S = Synapses(net,model="w : volt",on_pre="v_post+=w",method="euler")
    S.connect(p=p,condition="i!=j")
    return S.i[:],S.j[:]

def gen_unif_ds(Nsyns,dmin=.6,dmax=1.): return dmin+(dmax-dmin)*np.random.random(size=(1,Nsyns))

def gen_lognormal_gms(Nsyns,mean=0.,sigma=1.): return np.random.lognormal(mean=mean,sigma=sigma,size=(1,Nsyns))

def ConductAndReversPotsWithGapJunct(ConductOrig,ReversPotOrig,gLmin=2.,ggapHigh=1.3,full_dist=False,Ngaps=20,sigma=.4):
  NumNeurons = len(ConductOrig)
  ConductWithGapJunct = np.copy(ConductOrig)
  ReversPotWithGapJunctAux = ConductOrig*ReversPotOrig
  ReversPotWithGapJunct = np.zeros(np.shape(ReversPotOrig))
  synsgj = []
  ggs = []
  for presynGJ in range(NumNeurons):
    NumGapJuncts = 0
    NumIter = 0
    while NumGapJuncts <= int(Ngaps/2.) and NumIter<2*NumNeurons: # We start adding ggaps with Ngaps/2 partners (plus bidirectional, approx Ngaps partners) to each pre. We put a maximum in the iteration number, just in case the amount of possible gap junctions is less than Ngaps/2
      postsynGJ = np.random.randint(0,NumNeurons)
      if presynGJ==postsynGJ: continue # To avoid self-gap junctions
      if (presynGJ,postsynGJ) not in synsgj:
        if full_dist==True: ggap = [ sigma*truncnorm.rvs(a=.0,b=1.5/sigma) if np.random.rand()<.75 else ggapHigh ][0]
        if full_dist==False: ggap = sigma*truncnorm.rvs(a=.0,b=1.5/sigma)
        if ConductWithGapJunct[presynGJ]-ggap>=gLmin and ConductWithGapJunct[postsynGJ]-ggap>=gLmin:
          #Netpyne auto makes these junctions bidirectional so we don't want to read the direction in twice
          synsgj.append((presynGJ,postsynGJ)) # Bi-directional gap junctions
          #synsgj.append((postsynGJ,presynGJ)) # Bi-directional gap junctions
          ggs.append(ggap); #ggs.append(ggap) # Bi-directional gap junctions
          ConductWithGapJunct[presynGJ] -= ggap
          ConductWithGapJunct[postsynGJ] -= ggap # We reduce Leakage conductance
          ReversPotWithGapJunctAux[presynGJ] -= ggap*ReversPotOrig[postsynGJ]
          ReversPotWithGapJunctAux[postsynGJ] -= ggap*ReversPotOrig[presynGJ]
          NumGapJuncts+=1
      NumIter+=1

    ReversPotWithGapJunct[presynGJ] = ReversPotWithGapJunctAux[presynGJ]/ConductWithGapJunct[presynGJ]

  return ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj

def NetworkParams(NumNeurons=100,
FactorTau=1,
FactorKv3=1,
FactorKv7=1,
GapJunctProb=0.2,
ChemycalConnProb=0.36,
delaymin=.6,
delaymax=1.,
meangms=0.,
sigmagms=1.,
homogeneous=False,
NeuronModel=43,
randomseed = 7894 # FIX IT FOR REPRODUCIBILITY
):
  np.random.seed(randomseed)
  prefs.codegen.target = "numpy"
  seed(randomseed)
  try: 
    gLs = np.loadtxt("gLsLast.dat")
    ELs = np.loadtxt("ELsLast.dat")
    taums = FactorTau*np.loadtxt("taumsLast.dat")
    gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s = np.transpose(np.loadtxt( os.getcwd()+"ActiveParamsLast.dat" ))
  except:# Just in case you call the function form another directory
    print("Files not present in current directory. Searching in subdirectory") 
    gLs = np.loadtxt(os.getcwd()+"/Inet/gLsLast.dat")
    ELs = np.loadtxt(os.getcwd()+"/Inet/ELsLast.dat")
    taums = FactorTau*np.loadtxt(os.getcwd()+"/Inet/taumsLast.dat")
    gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s = np.transpose(np.loadtxt(os.getcwd()+"/Inet/ActiveParamsLast.dat" ))

  if homogeneous:
    gLs = [gLs[NeuronModel]]
    ELs = [ELs[NeuronModel]]
    taums = [taums[NeuronModel]]
    gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s = [gNas[NeuronModel]], [gKv3s[NeuronModel]], [gKv7s[NeuronModel]], [thm1s[NeuronModel]], [thh2s[NeuronModel]], [thn1s[NeuronModel]], [tha1s[NeuronModel]]

  NumModels = len(gLs) # Number of calibrated models
  
  if NumNeurons > NumModels:
    q, mod = divmod(NumNeurons, NumModels)
    gLs = np.concatenate((np.tile(gLs,q),gLs[:mod]))
    ELs = np.concatenate((np.tile(ELs,q),ELs[:mod]))
    taums = np.concatenate((np.tile(taums,q),taums[:mod]))
    gNas = np.concatenate((np.tile(gNas,q),gNas[:mod]))
    gKv3s = np.concatenate((np.tile(gKv3s,q),gKv3s[:mod]))
    gKv7s = np.concatenate((np.tile(gKv7s,q),gKv7s[:mod]))
    thm1s = np.concatenate((np.tile(thm1s,q),thm1s[:mod]))
    thh2s = np.concatenate((np.tile(thh2s,q),thh2s[:mod]))
    thn1s = np.concatenate((np.tile(thn1s,q),thn1s[:mod]))
    tha1s = np.concatenate((np.tile(tha1s,q),tha1s[:mod]))
  else:
    gLs = gLs[:NumNeurons]
    ELs = ELs[:NumNeurons]
    taums = taums[:NumNeurons]
    gNas = gNas[:NumNeurons]
    gKv3s = gKv3s[:NumNeurons]
    gKv7s = gKv7s[:NumNeurons]
    thm1s = thm1s[:NumNeurons]
    thh2s = thh2s[:NumNeurons]
    thn1s = thn1s[:NumNeurons]
    tha1s = tha1s[:NumNeurons]

  ENa=50; sigm1=4.; km2=.1; sigm2=13.; kh1=.012; kh2=.2; sigh1=-20.; sigh2=3.5; #Sodium current parameters
  EK=-90; sign1=12.; kn1 = FactorKv3*1.; kn2=FactorKv3*.001; sign2=-8.5; #Kv3 current parameters
  ka1 = FactorKv7*1; ka2=FactorKv7*.02; siga1=12.; siga2=-80.; #Kv7 current parameters (the rest of parameters are as in Kv3)
  SharedParams = [ENa, sigm1, km2, sigm2, kh1, kh2, sigh1, sigh2, EK, sign1, kn1, kn2, sign2, ka1, ka2, siga1, siga2]

  if GapJunctProb < 1e-3: 
    Ngaps=0
    ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj = np.copy(gLs), np.copy(ELs), [], []
  else:
    Ngaps = int(NumNeurons*GapJunctProb)
    ConductWithGapJunct, ReversPotWithGapJunct, ggs, synsgj = ConductAndReversPotsWithGapJunct(gLs, ELs,Ngaps=Ngaps)
  
  syns = gen_syns_samepop(NumNeurons,ChemycalConnProb)
  Nsyns = len(syns[0])  
  delays = gen_unif_ds(Nsyns,dmin=delaymin,dmax=delaymax)[0]
  gms = gen_lognormal_gms(Nsyns,mean=meangms,sigma=sigmagms)[0]

  syns = [(i,j) for (i,j) in zip(syns[0],syns[1])]

  CapsOrig =  taums*gLs/1e3
  CapsMod =  taums*ConductWithGapJunct/1e3
  #print(CapsOrig); quit()
  gms = [i*1e-3 for i in gms]
  delays = [i for i in delays]
  ggs = [i*1e-3 for i in ggs] 

  return gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, gNas, gKv3s, gKv7s, thm1s, thh2s, thn1s, tha1s, SharedParams, syns, delays, gms, synsgj, ggs

### USAGE
#NumNetw = 2
#for k in range(NumNetw):
#  gLs, ELs, CapsOrig, ConductWithGapJunct, ReversPotWithGapJunct, CapsMod, syns, synsgj, SharedParams = NetParams(NumNeurons=100, FactorTau=1, FactorKv3=1, FactorKv7=1, GapJunctProb=0.1, ChemycalConnProb=0.36, delaymin=.6, delaymax=1., meangms=0., sigmagms=1.)

#  print(gLs[:2], ELs[:2], CapsOrig[:2], ConductWithGapJunct[:2], ReversPotWithGapJunct[:2], CapsMod[:2])








