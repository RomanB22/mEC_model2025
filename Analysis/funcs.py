# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import cwt, morlet2, butter, sosfiltfilt, welch
import matplotlib.pyplot as plt

def calculate_psd(signal, dt_ms):
    """Calculates the Power Spectral Density (PSD) using Welch's method."""
    dt_sec = dt_ms * 1e-3  # Convert milliseconds to seconds
    fs = 1 / dt_sec  # Sampling frequency

    # Using Welch's method for PSD estimation (more robust than direct FFT)
    frequencies, psd = welch(signal, fs, nperseg=4096)  # Adjust nperseg as needed

    return frequencies, psd

def find_peak_frequency(frequencies, psd):
    """Finds the frequency at the peak of the PSD."""
    peak_index = np.argmax(psd)
    peak_frequency = frequencies[peak_index]
    return peak_frequency

def compPWT(signal,fs=np.arange(50.,450.1,3.),dt=1.e-4,word=5,tmin=.0):
  """
  Computes the Wavelet Power Spectra (modulus of complex wavelet transform) for a set of signals
  Input:
   signal (to be analyzed)
   fs (frequencies to compute cwt)
   dt
   word (wavelet order, default=5)
   nlevels (number of levels at contour plot)
   tmin (plot wavelet spectra from tmin to end of signal. in seconds)
  """

  ktmin = int(tmin/dt);

  # Perform Wavelet Analysis
  samp_rate = 1./dt; # Sampling rate is inverse of time step dt
  ss = word*samp_rate/(2.*np.pi*fs); # Compute corresponding wavelet scalings or widths (s).
  pwt = np.power(np.abs(cwt(signal[ktmin:]-np.mean(signal[ktmin:]),morlet2,widths=ss)),2.) # Compute wavelet power, i.e. squared modulus, using a Morlet wavelet

  return pwt

def compPhaseWT(signal,fs=np.arange(50.,450.1,3.),dt=1.e-4,word=5,tmin=.0):
  """
  Computes the Wavelet Power Phase (modulus of complex wavelet transform) for a set of signals
  Input:
   signal (to be analyzed)
   fs (frequencies to compute cwt)
   dt
   word (wavelet order, default=5)
   nlevels (number of levels at contour plot)
   tmin (plot wavelet spectra from tmin to end of signal. in seconds)
  """

  ktmin = int(tmin/dt);

  # Perform Wavelet Analysis
  samp_rate = 1./dt; # Sampling rate is inverse of time step dt
  ss = word*samp_rate/(2.*np.pi*fs); # Compute corresponding wavelet scalings or widths (s).
  #pwt = np.angle(cwt(signal[ktmin:]-np.mean(signal[ktmin:]),morlet2,widths=ss)) # Compute wavelet phase, using a Morlet wavelet
  pwt = np.angle(cwt(signal[ktmin:]-np.mean(signal[ktmin:]),morlet2,widths=ss))
  return pwt

def plotPow(cwtPow,max_power=None,fs=np.arange(50.,450.1,3.),tmin=.0,dt=.1,levels=[0,1],nlevels=30,Ncycs=1,flims=[],xlims=[],colorbar=True,cmap="hot",thetaPeriod=125.):
    """
    plot the contour plot for the wavelet scalogram for one single simulation
    """
    ktmin = int(np.round(tmin/dt));
    Nt_pr_cyc = int(np.round(thetaPeriod/dt));
    if max_power==None: max_power = np.max(np.max(cwtPow));
    #levels = np.linspace(0.,max_power+1.,nlevels);
    print('Printing length of inputs to contourf',np.shape(np.linspace(-np.pi,-np.pi+2.*np.pi*Ncycs,int(np.round(Ncycs*Nt_pr_cyc)))), np.shape(fs), np.shape(np.abs(cwtPow[:,ktmin:])), len(levels))
    tmp = plt.contourf(np.linspace(-np.pi,-np.pi+2.*np.pi*Ncycs,int(np.round(Ncycs*Nt_pr_cyc))),fs,np.abs(cwtPow[:,ktmin:]),levels=levels,cmap=cmap,rasterized=True);

    if colorbar==True: tmp = plt.colorbar(tmp);

    """if len(flims)==0: tmp = ax.set_ylim(fs[0],fs[-1]);
    else: tmp = ax.set_ylim(flims[0],flims[1]);

    if len(xlims)==0: tmp = ax.set_xlim(0,dt*np.shape(cwtPow)[1]);
    else: tmp = ax.set_xlim(xlims[0],xlims[1]);"""

    if len(flims)==0: tmp = plt.ylim(fs[0],fs[-1]);
    else: tmp = plt.ylim(flims[0],flims[1]);

    if len(xlims)==0: tmp = plt.xlim(-np.pi,-np.pi+2.*np.pi*Ncycs);
    else: tmp = plt.xlim(xlims[0],xlims[1]);

def plotPhase(cwtPow,max_power=None,fs=np.arange(50.,450.1,3.),tmin=.0,dt=.1,levels=[0,1],nlevels=30,Ncycs=1,flims=[],xlims=[],colorbar=True,cmap="hot",thetaPeriod=125.):
    """
    plot the contour plot for the wavelet scalogram for one single simulation
    """
    ktmin = int(np.round(tmin/dt));
    Nt_pr_cyc = int(np.round(thetaPeriod/dt));
    #if max_power==None: max_power = np.max(np.max(cwtPow));
    levels = 30#np.linspace(0.,max_power+1.,nlevels);
    tmp = plt.contourf(np.linspace(-np.pi,-np.pi+2.*np.pi*Ncycs,int(np.round(Ncycs*Nt_pr_cyc))),fs,cwtPow[:,ktmin:],levels=levels,cmap=cmap,rasterized=True);

    if colorbar==True: tmp = plt.colorbar();

    if len(flims)==0: tmp = plt.ylim(fs[0],fs[-1]);
    else: tmp = plt.ylim(flims[0],flims[1]);

    if len(xlims)==0: tmp = plt.xlim(-np.pi,-np.pi+2.*np.pi*Ncycs);
    else: tmp = plt.xlim(xlims[0],xlims[1]);

def plotMaskedPhase(cwtPow,cwtPhase,max_power=None,fs=np.arange(50.,450.1,3.),tmin=.0,dt=.1,levels=[0,1],nlevels=30,Ncycs=1,flims=[],xlims=[],colorbar=True,cmap="hot",thetaPeriod=125.):
    """
    plot the contour plot for the phase with power 2 times standard deviation above mean
    """
    ktmin = int(np.round(tmin/dt));
    Nt_pr_cyc = int(np.round(thetaPeriod/dt));
    Auxcwt = cwtPhase[:,ktmin:]

    #if max_power==None: max_power = np.max(np.max(cwtPow));
    levels = 30#np.linspace(0.,max_power+1.,nlevels);
    mask2D = np.ma.MaskedArray(cwtPhase, mask=cwtPow > np.mean(cwtPow) + 2*np.std(cwtPow))
    
    tmp = plt.contourf(np.linspace(-np.pi,-np.pi+2.*np.pi*Ncycs,int(np.round(Ncycs*Nt_pr_cyc))),fs,Auxcwt *np.ma.getmaskarray(mask2D),levels=levels,cmap=cmap,rasterized=True);

    if colorbar==True: tmp = plt.colorbar();

    if len(flims)==0: tmp = plt.ylim(fs[0],fs[-1]);
    else: tmp = plt.ylim(flims[0],flims[1]);

    if len(xlims)==0: tmp = plt.xlim(-np.pi,-np.pi+2.*np.pi*Ncycs);
    else: tmp = plt.xlim(xlims[0],xlims[1]);


