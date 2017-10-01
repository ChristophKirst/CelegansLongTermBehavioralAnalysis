# -*- coding: utf-8 -*-
"""
Frequency analysis

@author: ckirst
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


import experiment as exp;
import plot as fplt;


strain = 'N2';
feat = ['roam', 'speed'];

data = exp.load_data(strain);

iworm = 0;
ntimes = data.total_time[iworm];
dt = 1.0 / 3.0; # 3 Hz time resolution

for fi,fn in enumerate(feat):
  d = getattr(data, fn);
  
  #rf = np.fft.fft(r[iworm, :ntimes]);
  #rf = 2.0/ntimes * np.abs(rf[:ntimes//2]);
  #ff = np.linspace(0.0, 1.0/(2.0*dt), ntimes//2)

  ### FFT
  ff, pfsd_r, = signal.periodogram(d[iworm, :ntimes], 1/dt);

  nmean = 50;
  pfsd_r_m = np.convolve(pfsd_r, np.ones(nmean)/nmean, mode='same');
  
  plt.figure(1*fi); plt.clf();
  plt.semilogy(ff, pfsd_r)
  plt.semilogy(ff, pfsd_r_m)
  plt.tile('power %s %s %d' % (strain, fn, iworm));
  
  
  ### Spectogram
  f, t, Sxx = signal.spectrogram(d[iworm,:ntimes], 1/dt, nperseg = 2**14, noverlap = 2**14-2**8);
  
  plt.figure(2*fi); plt.clf();
  ax = plt.subplot(2,1,1);
  nfmax = 50;
  plt.pcolormesh(t, f[:nfmax], Sxx[:nfmax,:])
  plt.ylabel('Frequency [Hz]')
  plt.xlabel('Time [sec]')
  plt.tile('spectogram %s %s %d' % (strain, fn, iworm));
  plt.subplot(2,1,2, sharex = ax);
  tt = np.linspace(0, t[-1], ntimes);
  plt.plot(tt, d[iworm,:ntimes]);


  plt.figure(3*fi); plt.clf();
  ax = plt.subplot(2,1,1);
  smax = 0.05;
  nfmax = -1;
  sxx2 = Sxx.copy();
  sxx2[sxx2 > smax] = smax;
  plt.pcolormesh(t, f[:nfmax], sxx2[:nfmax,:])
  plt.ylabel('Frequency [Hz]')
  plt.xlabel('Time [sec]')
  plt.show()
  plt.subplot(2,1,2, sharex = ax);
  tt = np.linspace(0, t[-1], ntimes);
  plt.plot(tt, d[iworm,:ntimes]);
  
  fplt.plot_array(sxx2, title = 'spectorgram %s %s %d' % (strain, fn, iworm))
  
  plt.figure(4*fi); plt.clf();
  plt.plot(np.sum(sxx2, axis = 0))
  plt.title('sum spectorgram %s %s %d' % (strain, fn, iworm))
