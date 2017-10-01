# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:21:32 2016

@author: ckirst
"""

import numpy as np
#import matplotlib.pyplot as plt
#import scipy.signal as sig

import experiment as exp;
import analysis as alys;
import plot as fplt

### Load Data

strains = ['cat-2', 'tdc-1', 'daf-7', 'tph-1']
strain = 'N2';

for strain in strains:
  feat = 'roam';
  
  save_fig = True;
  
  data = exp.load_data(strain);
  
  ### Averaging at differnet scales
  
  def kernel_flat(s):
    return 1.0/s * np.ones(s);
    
  def kernel_gauss(s, nsigma = 3):
    x = np.linspace(-nsigma, nsigma, s);
    k = np.exp(-x*x/2);
    return k / k.sum();
  
  def average(data, kernel, scales):
    res = np.zeros((data.shape[0], len(scales), data.shape[1]));
    for si,s in enumerate(scales):
      kernel_data = kernel(s);
      for i in range(data.shape[0]):
        res[i,si] = np.convolve(data[i], kernel_data, mode = 'same');
    return res;
  
  ### Run parallel
  
  # write data to npy for memmap use
  d = getattr(data, feat);
  shape = d.shape;
  f = exp.save(d, '%s_%s.npy' % (strain, feat));
  
  ### Create numpy file for results
  
  scales = [2**i+1 for i in range(1,17)];
  
  res_shape = (shape[0], len(scales), shape[1]);
  exp.memmap('%s_%s_scales_mean.npy' % (strain, feat), dtype = 'float32', mode = 'w+', shape = res_shape)
  
  
  def run(i):
    print 'processing worm %d' % i
    fdat = exp.load('%s_%s.npy' % (strain, feat), memmap = 'r');
    signal = fdat[[i],:];
    res = average(signal, kernel_flat, scales);
    #res = average(signal, kernel_gauss, scales);
    fres = exp.load('%s_%s_scales_mean.npy' % (strain, feat), memmap = 'r+');
    fres[i,:,:] = res[0];
  
  from multiprocessing import Pool #, cpu_count;
  
  #pool = Pool(processes = cpu_count()-6);
  pool = Pool(processes = 12);
  
  pool.map(run, range(data.nworms));
  
  
  ### Time normalize data
  fres = exp.load('%s_%s_scales_mean.npy' % (strain, feat), memmap = 'r');
  
  sbins = exp.stage_bins(data, nbins = 8192);
  res_shape = (fres.shape[0], fres.shape[1], sbins[0].shape[1]);
  
  exp.memmap('%s_%s_time_normalized_scales_mean.npy' % (strain, feat), shape = res_shape, mode = 'w+');
  
  def norm(i):
    print 'scale %d' % i;
    fres    = exp.load('%s_%s_scales_mean.npy' % (strain, feat), memmap = 'r');
    fres_tn = exp.load('%s_%s_time_normalized_scales_mean.npy' % (strain, feat), memmap = 'r+');
    fres_tn[:,i,:] = exp.bin_data(fres[:,i,:], sbins);
  
  pool = Pool(processes = 12);
  
  pool.map(norm, range(fres.shape[1]));
  
  
  # plot
  res_tn = exp.load('%s_%s_time_normalized_scales_mean.npy' % (strain, feat), memmap = 'r');
  
  dorder = exp.load('%s_%s_order.npy' % (strain, feat))
  
  datplt = alys.scales_to_array(res_tn, worms_first=False, order = dorder);
  fig = fplt.plot_array(datplt, title = '%s %s time normalized scales mean' % (strain, feat));
  fplt.savefig(fig,exp.figname('%s_%s_time_normalized_scales_mean.png'% (strain, feat)), width = 2500)
  
  datplt = alys.scales_to_array(res_tn, worms_first=True, order = dorder);
  fig = fplt.plot_array(datplt, title = '%s %s time normalized scales mean' % (strain, feat));
  fplt.savefig(fig,exp.figname('%s_%s_time_normalized_scales_mean_2.png'% (strain, feat)), width = 2500)
  
  
  #### Test on single worm
  # 
  #from timer import timeit;
  #     
  #scales = [2**i+1 for i in range(1,14)]
  #print scales
  #
  #speed_test = data.speed[:10, :];
  #
  #@timeit
  #def run_test():
  #  return average(speed_test, kernel_flat, scales);
  #
  #av = run_test();
  #
  #plt.figure(1); plt.clf();
  #nplt = av.shape[0]
  #for i in range(nplt):
  #  plt.subplot(nplt,1,i+1);
  #  plt.imshow(av[i], aspect='auto', cmap='jet', vmin = 0)
  #plt.tight_layout();
  #
  #
  #
  #### Test on roaming and dwelling
  #
  #scales = [2**i+1 for i in range(1,14)]
  #print scales
  #
  #test = data.roam[:10, :];
  #
  #@timeit
  #def run_test():
  #  return average(test, kernel_flat, scales);
  #
  #av = run_test();
  #
  #plt.figure(2); plt.clf();
  #nplt = av.shape[0]
  #for i in range(nplt):
  #  plt.subplot(nplt,1,i+1);
  #  d = av[i,:];
  #  for i in range(d.shape[0]):
  #    #d[i] /= (np.median(d[i])+1.0);
  #    d[i] /= d[i].max();
  #    d[i][d[i] > 0.9] = 0.9;
  #  vmax = d.max();
  #  for k in range(-200,200,1):
  #      d[:,data.stage_switch[i] + k] = vmax;
  #  plt.imshow(d, aspect='auto', cmap='jet', vmin = 0)
  #plt.tight_layout();
