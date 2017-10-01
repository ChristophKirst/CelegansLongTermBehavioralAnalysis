# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 15:22:21 2016

@author: ckirst
"""

import os

dir_other = '/home/ckirst/Science/Projects/CElegansBehaviour/Experiment/DwellingRoaming/Scripts/Other/'
dir_base = '/home/ckirst/Science/Projects/CElegansBehaviour/Experiment/DwellingRoaming'

os.chdir(dir_other);

import pykov as pk

import numpy as np


t = np.random.rand(100) > 0.5;
t = np.array(t, dtype = int);

p,c = pk.maximum_likelihood_probabilities(t)


os.chdir(dir_base);

import experiment as exp;


strain = 'N2';
feat = 'roam';

data = exp.load_data(strain);
d = getattr(data, feat)


sbins = exp.stage_bins(data, nbins = 2**4);


c = np.zeros((data.nworms, sbins[0].shape[1], 2));
p = np.zeros((data.nworms, sbins[0].shape[1], 2));
for wid in range(data.nworms):
  print 'worm %d' % wid
  st = sbins[0][wid];
  ed = sbins[1][wid];
  for i,se in enumerate(zip(st,ed)):
    s,e = se;
    pp,cc = pk.maximum_likelihood_probabilities(d[wid,s:e]);
    c[wid, i] = [cc[(0,1)], cc[(1,0)]];
    p[wid, i] = [pp[0], pp[1]];

    
### plot the result

import matplotlib.pyplot as plt
plt.figure(1); plt.clf();
wid = 0;
plt.subplot(2,1,1);
plt.plot(c[wid,:,0]);
plt.plot(c[wid,:,1]);
plt.subplot(2,1,2);
plt.plot(p[wid,:,0]);
plt.plot(p[wid,:,1]);


import plot as fplt




fplt.plot_array(p[:,:,1], title = '%s %s probability' % (strain, feat))
#fplt.plot_array(p[:,:,1])

cc = c.copy();
cc[cc > 0.5] = 0.5;

fplt.plot_array(cc[:,:,0], title = '%s %s 0->1' % (strain, feat))
fplt.plot_array(cc[:,:,1], title = '%s %s 1->0' % (strain, feat))  


