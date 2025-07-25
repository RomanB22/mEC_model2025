import numpy as np
from scipy.signal import cwt, morlet2, butter, sosfiltfilt
import matplotlib.pyplot as plt

def compPWT(signal,dt,fs=np.arange(50.,450.1,3.),word=5):
  """
  Computes the Wavelet Power Spectra (modulus of complex wavelet transform) and Phase for a set of signals
  Input:
   signal (to be analyzed)
   fs (frequencies to compute cwt)
   sampling rate
   word (wavelet order, default=5)
  """
  samp_rate=1./dt
  # Perform Wavelet Analysis
  ss = word*samp_rate/(2.*np.pi*fs); # Compute corresponding wavelet scalings or widths (s).
  prwt = np.power(np.abs(cwt(signal-np.mean(signal),morlet2,widths=ss)),2.) # Compute wavelet power, i.e. squared modulus, using a Morlet wavelet

  phwt = np.angle(cwt(signal-np.mean(signal),morlet2,widths=ss))  # Compute wavelet phase, using a Morlet wavelet
  return prwt, phwt

def plotPow(cwtPow,max_power=None,fs=np.arange(50.,450.1,3.),dt=.1,levels=None,nlevels=30,flims=[],xlims=[],colorbar=True,cmap="hot",ax=None,fig=None):
    """
    plot the contour plot for the wavelet scalogram for one single simulation
    """

    if max_power==None and levels is None: max_power = np.max(np.max(cwtPow)); levels = np.linspace(0., max(max_power,0.1), nlevels, endpoint=True)
    tmp = ax.contourf(np.linspace(0.,dt*np.shape(cwtPow)[1],np.shape(cwtPow)[1],endpoint=True),fs,np.abs(cwtPow),levels=levels,cmap=cmap,rasterized=True);

    if colorbar==True: tmp = fig.colorbar(tmp, ax=ax);

    if len(flims)==0: tmp = ax.set_ylim(fs[0],fs[-1]);
    else: tmp = ax.set_ylim(flims[0],flims[1]);

    if len(xlims)==0: tmp = ax.set_xlim(0,dt*np.shape(cwtPow)[1]);
    else: tmp = ax.set_xlim(xlims[0],xlims[1]);

def plotMaskedPhase(cwtPow,cwtPhase,max_power=None,fs=np.arange(50.,450.1,3.),dt=.1,levels=[0,1],nlevels=30,flims=[],xlims=[],colorbar=True,cmap="hot",ax=None,fig=None):
    """
    plot the contour plot for the phase with power 2 times standard deviation above mean
    """
    Auxcwt = cwtPhase

    mask2D = np.ma.MaskedArray(cwtPhase, mask=cwtPow > np.mean(cwtPow) + 2*np.std(cwtPow))
    
    tmp = ax.contourf(np.linspace(0,np.shape(cwtPow)[1]*dt,np.shape(cwtPow)[1],endpoint=True),fs,Auxcwt *np.ma.getmaskarray(mask2D),levels=nlevels,cmap=cmap,rasterized=True);

    if colorbar==True: tmp = fig.colorbar(tmp);

    if len(flims)==0: tmp = ax.set_ylim(fs[0],fs[-1]);
    else: tmp = ax.set_ylim(flims[0],flims[1]);

    if len(xlims)==0: ax.set_xlim(0,dt*np.shape(cwtPhase)[1]);
    else: tmp = ax.set_xlim(xlims[0],xlims[1]);

