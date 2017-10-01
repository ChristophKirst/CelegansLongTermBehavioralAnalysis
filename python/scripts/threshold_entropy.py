# -*- coding: utf-8 -*-
"""
Stereotypy analysis 

How stereotypic is the behavior at different time scales

The analysis works as follows:
  * threshold the time normalized and scale averaged data to get active inactive regions on different scales
  * 1-variance at each time bin, time scale and threshold then gives a measure of 
"""

### threshold + variance to identify consistent activity perids art diferrent time scales

import numpy as np
import matplotlib.pyplot as plt

import scipy.signal as sig

import experiment as exp;
import analysis as alys
import plot as fplt


strain = 'tdc-1';
feat = 'roam';

data = exp.load('%s_%s_time_normalized_scales_mean.npy' % (strain, feat))
nworms, nscales, ntimes = data.shape

data_ord = exp.load('%s_%s_order.npy' % (strain, feat));


plt.figure(70); plt.clf();
plt.hist(data[0].T, bins = 50, histtype = 'bar')
plt.title('%s %s time normalized scales mean distribution' % (strain, feat))

### Threshold use multiple thresholds per scale

thres = np.linspace(0,1,20)[1:-1];
nthres = len(thres);

data_th = np.zeros((nworms, nthres, nscales, ntimes));
for i,th in enumerate(thres):
  data_th[:,i,:,:] =  data > th;


#for ith in range(nthres):
#  datp = alys.scales_to_array(data_th[:,ith,:,:], worms_first = True);
#  fig = fplt.plot_array(datp, title = 'thresholded th = %f' % thres[ith]);
#  fplt.savefig(fig, exp.figname('%s_%s_time_normalized_scales_mean_thresholded_th=%f.png' % (strain, feat, thres[ith])), width = 2500)
  

### Std along worm axis 
# Note: for 0,1 values std = mean (1-mean)

data_th_mean = np.mean(data_th, axis  = 0);
data_th_std  = np.std( data_th, axis  = 0);

dorder = exp.load('%s_%s_order.npy' % (strain, feat))

datp = alys.scales_to_array(data_th_std, worms_first = True, );
fig = fplt.plot_array(datp, title = '%s %s thresholded std' % (strain, feat));
fplt.savefig(fig, exp.figname('%s_%s_time_normalized_scales_mean_thresholded_std.png' % (strain, feat)), width = 2500)

#plt.figure(7); plt.clf();
#plt.imshow(datp[::-1,:], aspect = 'auto', extent = [0, ntimes, 0, datp.shape[1]])
#plt.tight_layout();

datp = alys.scales_to_array(data_th_std, worms_first = False);
fplt.plot_array(datp, title = '%s %s thresholded std' % (strain, feat));

#fig = plt.figure(8); plt.clf();
#plt.imshow(datp[::-1,:], aspect = 'auto', extent = [0, ntimes, 0, datp.shape[1]])
#plt.tight_layout();



###

#datp = exp.scales_to_array(data_th_mean, worms_first = True);
#fplt.plot_array(datp, title = 'thresholded mean');
#
#plt.figure(70); plt.clf();
#plt.imshow(datp[::-1,:], aspect = 'auto', extent = [0, ntimes, 0, datp.shape[1]])
#plt.tight_layout();
#
#
#datp = exp.scales_to_array(data_th_mean, worms_first = False);
#fplt.plot_array(datp, title = 'thresholded mean');
#
#plt.figure(80); plt.clf();
#plt.imshow(datp[::-1,:], aspect = 'auto', extent = [0, ntimes, 0, datp.shape[1]])
#plt.tight_layout();



#plt.figure(300); plt.clf();
#def red(x):
#  return [x, 0.0, 0.0, 1.0];
#for i in range(nscales):
#  plt.plot(data_th_std[0,i,:], linewidth = 2, color = red(i*1.0/(nscales-1)) )
  
ith = 0;
fig = plt.figure(301); plt.clf();
for s in range(nscales):
  plt.subplot(4,4,s+1)
  plt.plot((0.5-data_th_std[ith,s, :])*2, linewidth = 1, color = 'r' )
  plt.plot(data_th_mean[ith,s,:], linewidth= 1, color = 'k' )
  plt.title('%s %s scale %d' % (strain, feat, (2**(i+1))));
plt.tight_layout();

fig.savefig(exp.figname('%s_%s_time_normalized_scales_mean_thresholded_th=%f_stereotypy.png' % (strain, feat, thres[ith])));


ith = 1;
fig = plt.figure(302); plt.clf();
for i in range(nscales):
  plt.subplot(4,4,i+1)
  plt.plot((0.5-data_th_std[0,ith,:])*2, linewidth = 1, color = 'r' )
  plt.plot(data_th_mean[0,ith,:], linewidth= 1, color = 'k' )
  plt.title('%s %s scale %d' % (strain, feat, (2**(i+1))));
plt.tight_layout();

fig.savefig(exp.figname('%s_%s_time_normalized_scales_mean_thresholded_th=%f_stereotypy.png' % (strain, feat, thres[ith])));