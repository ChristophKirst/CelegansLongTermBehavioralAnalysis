# -*- coding: utf-8 -*-
"""
Data Preparation and Basic Analysis for the Roaming and Dwelling data sets

@author: ckirst
"""
import os
import numpy as np
import matplotlib.pyplot as plt

import experiment as exp

#
#strains = exp.strains;
#strains = ['N2']
#for strain in strains:
#  exp.convert_data(strain)


strains = exp.strains;
strains = [ 'tdc-1', 'daf-7', 'tph-1']
strains = ['N2'];

strain = 'N2';
feat = 'roam';


save_fig = exp.figname('2016_12_20')
save_fig = None;
  
for strain in strains:
  ### Load worms
  
  data = exp.load_data(strain);
  nworms = data.nworms;
  
  ### Stage time distribution
  
  ssort = np.argsort(data.total_time);
  
  fig = plt.figure(1); plt.clf();
  plt.subplot(3,2,1)
  plt.plot(np.array(data.total_time)[ssort]);
  #for s in range(data.nstages-1):
  #  o = np.argsort(data.stage_durations[:,s]);
  #  plt.plot(np.array(data.total_time)[o]);
  plt.title('%s %s - total time' % (strain, feat));
  
  for s in range(data.nstages):
    plt.subplot(3,2,s+2);
    plt.plot(data.stage_durations[ssort,s]);
    #for si in range(data.nstages-1):
    #  plt.plot(data.stage_durations[np.argsort(data.stage_durations[:,si]), s]);
    #plt.title('mean in L%d' % (s+1));
  plt.tight_layout();
  
  if save_fig is not None:
    fig.savefig(os.path.join(save_fig, '%s_%s.png' % (strain, 'stage_durations')))
  
  
  ### Mean feature accross worms and stages
  
  d = getattr(data, feat);
  
  #average+order
  dm = np.sum(d, axis = 1) * 1.0 / data.total_time;
  
  dorder = np.argsort(dm);
  exp.save(dorder, '%s_%s_order.npy' % (strain, feat));
  
  #stage resolved averages
  sbins = exp.stage_bins(data, nbins = 1)
  dms = exp.bin_data(d, sbins)
  dms_shuffle = dms.copy();
  for i in range(dms.shape[1]):
    dms_shuffle[:,i] = np.random.permutation(dms_shuffle[:,i]);
  dms_data = [dms, dms_shuffle];
  
  
  
  fig = plt.figure(20+i); plt.clf();
  plt.subplot(3,2,1)
  plt.plot(dm[dorder]);
  for s in range(data.nstages):
    plt.plot(dm[np.argsort(dms[:,s])]);
  plt.title('%s %s - total mean' % (strain, feat));
  
  for s in range(data.nstages):
    plt.subplot(3,2,s+2);
    plt.plot(dms[dorder,s]);
    for si in range(data.nstages):
      plt.plot(dms[np.argsort(dms[:,si]), s]);
    plt.title('mean in L%d' % (s+1));
  plt.tight_layout();
  
  if save_fig is not None:
    fig.savefig(os.path.join(save_fig,'%s_%s_basic_analysis.png' % (strain, feat)))
  
  
  
  ### PCA on feature means
  
  from matplotlib.mlab import PCA
  from mpl_toolkits.mplot3d import Axes3D
  
  for i in range(len(dms_data)):
    ddms = dms_data[i];
    results = PCA(ddms);
    pcs = results.Y;
   
    fig = plt.figure(101+i); plt.clf();
    plt.subplot(1,3,1);
    plt.imshow(pcs[dorder], interpolation = 'none', aspect = 'auto', cmap = 'viridis')
    plt.title('pca components');
    plt.subplot(2,3,2);
    plt.imshow(results.Wt, cmap = 'magma', interpolation = 'none' );
    plt.colorbar()
    plt.title('pca vectors')
    ax = fig.add_subplot(2,3,5, projection = '3d');
    ax.scatter(pcs[:,0], pcs[:,1], pcs[:,2], 'bo');
    plt.xlabel('PCA1'); plt.ylabel('PCA2');
    ax.set_zlabel('PCA3');
    plt.subplot(2,3,3);
    plt.plot(results.mu)
    plt.title('mean');
    plt.subplot(2,3,6);
    plt.plot(np.cumsum(results.fracs), 'r')
    plt.title('variance explained')
    plt.ylim(0,1)
    plt.tight_layout();
    
    #if save_fig:
    #  fig.savefig(exp.figname('%s_%s_basic_analysis_pca.png' % (strain, feat)))
    
    # t - SNE on romaing fraction
    import plot
    reload(plot)
  
    plt.figure(11+i); plt.clf();
    Y = plot.plot_tsne(ddms)
    plt.colorbar()
  
  
    # find the outliers 
    ctr = np.mean(Y, axis = 0);
    ddist = np.linalg.norm(Y - ctr, axis = 1);
  
  
  for i in range(len(dms_data)):
    ddms = dms_data[i];
    # find outliers w.r.t mean
    ctr = np.mean(ddms, axis = 0);
    ddist = np.linalg.norm(ddms - ctr, axis = 1);
    
  
    plt.figure(60+i); plt.clf();
    plt.hist(ddist, 32)
  
  
  
  
  #plt.plot(ddist);
  
  
  #plt.scatter(Y[51:52,0], Y[51:52,1]);
  
  
  #outliers are 51 and 20 
  
  plt.figure(13); plt.clf();
  plt.subplot(1,3,1);
  plt.plot(dms[[20,51],:].T)
  plt.subplot(1,3,2);
  plt.imshow(dms);
  plt.subplot(1,3,3);
  plt.imshow(dms[dorder,:]) 
  
  ### PCA on switching and fraction
  
  plt.figure(200);
  plt.clf();
  
  
  fd = np.zeros((nworms,5));
  sn = np.zeros((nworms,5));
  for s in range(5):
    for i in dorder:
      dd = d[i];
      idx = data.stage[i] == s+1;
      fd[i,s] = dd[idx].sum() * 1.0 / len(idx);
      sn[i,s] = np.sum(np.diff(dd[idx])==1) * 1.0 / len(idx);
    
    plt.subplot(2,5,s+1);
    plt.plot(fd[:,s]);
    plt.title('fraction in L%d' % (s+1));
    
    plt.subplot(2,5,s+1+5);
    plt.plot(sn[:,s]);
    plt.title('fraction of switches in L%d' % (s+1));
  
  
  # pca
  
  a = np.hstack((fd, sn))
  
  results = PCA(a);
  pcs = results.Y;
   
  fig = plt.figure(101); plt.clf();
  plt.subplot(1,3,1);
  plt.imshow(pcs[dorder], interpolation = 'none', aspect = 'auto', cmap = 'viridis')
  plt.title('pca components');
  plt.subplot(2,3,2);
  plt.imshow(results.Wt, cmap = 'magma', interpolation = 'none' );
  plt.title('pca vectors')
  ax = fig.add_subplot(2,3,5, projection = '3d');
  ax.scatter(pcs[:,0], pcs[:,1], pcs[:,2], 'bo');
  plt.xlabel('PCA1'); plt.ylabel('PCA2');
  ax.set_zlabel('PCA3');
  plt.subplot(2,3,3);
  plt.plot(results.mu)
  plt.title('mean');
  plt.subplot(2,3,6);
  plt.plot(np.cumsum(results.fracs), 'r')
  plt.title('variance explained')
  plt.tight_layout();
  
  if save_fig:
    fig.savefig(os.path.join(save_fig, '%s_%s_switching_basic_analysis_pca.png' % (strain, feat)))
  
  
