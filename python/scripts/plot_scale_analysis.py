# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 12:33:16 2016

@author: ckirst
"""
import numpy as np
import experiment as exp
import analysis as alys
import plot as fplt


feat = 'roam';

nmax = 25;

for strain in exp.strains:
  strain = 'N2';
  #strain = 'tdc-1'
  res_tn = exp.load('%s_%s_time_normalized_scales_mean.npy' % (strain, feat), memmap = 'r');  
  dat = res_tn.copy();
  dorder = exp.load('%s_%s_order.npy' % (strain, feat))  
  #datplt = alys.scales_to_array(res_tn, worms_first=False, order = dorder);
  #n = min(nmax, res_tn.shape[0]);
  #idx = range(n);
  idx = dorder[25:50];
  for s in range(res_tn.shape[1]):
    dat[:,s,:] = (1.0 /  np.max(dat[:,s,:], axis = 1) *  dat[:,s,:].T).T ;
  datplt = alys.scales_to_array(dat[idx, 6:], worms_first=False)
  fig = fplt.plot_array(datplt, title = '%s %s time normalized scales mean' % (strain, feat));
  
  
  strain = 'N2';
  #strain = 'tdc-1'
  res_tn = exp.load('%s_%s_time_normalized_scales_mean.npy' % (strain, feat), memmap = 'r');  
  dat = res_tn.copy();
  dorder = exp.load('%s_%s_order.npy' % (strain, feat))  
  #datplt = alys.scales_to_array(res_tn, worms_first=False, order = dorder);
  #n = min(nmax, res_tn.shape[0]);
  #idx = range(n);
  idx = dorder[25:50];
  idx = range(res_tn.shape[0]);
  for s in range(res_tn.shape[1]):
    dat[:,s,:] = (1.0 /  np.max(dat[:,s,:], axis = 1) *  dat[:,s,:].T).T ;
  dat = dat > 0;
  datplt = alys.scales_to_array(dat[idx, 6:], worms_first=False)
  fig = fplt.plot_array(datplt, title = '%s %s time normalized scales mean' % (strain, feat));