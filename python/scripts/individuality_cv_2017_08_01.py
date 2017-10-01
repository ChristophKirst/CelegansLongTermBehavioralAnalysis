# -*- coding: utf-8 -*-
"""
Distribution analysis to identify individuality
"""


import os
import numpy as  np
import matplotlib.pyplot as plt;

import scipy.stats as stats

import experiment as exp
import analysis as als
import plot as fplt;


#%% Parameter

strains = exp.strains;

#strains = ['N2', 'ser-1', 'ser-4', 'ser-5', 'ser-7', 'mod-1', 'mod-5'];
strains = ['N2', 'tph-1', 'cat-2', 'tdc-1', 'npr-1', 'daf-7'];

strains_label = '_'.join(strains);



feat = 'speed';

nstrains = len(strains);
nstages = 5;

dbins = 10;

sbins = 37;
sbins = 20;
sbins = 75;
sbins = 10;
#sbins = 5;
sbins = 1;

tbins = nstages * sbins;


fig_dir = '/data/Science/Projects/CElegans/Analysis/Roaming/Figures/2017_08_01'
#fig_dir = None;

d = als.jensen_shannon_divergence;

def d(p,q):
  dd = als.jensen_shannon_divergence(p,q);
  if np.sum(np.linspace(0,1, len(p)) * p) > np.sum(np.linspace(0,1, len(q)) * q):
    return dd;
  else:
    return -dd;

#%% Load / Prepare data

data = {}; dat_bin = {}; dat_mean = {}; dat_var = {};
dat = {}; stage_bins = {}; dat_bin_s = {};

for strain in strains:
  print 'processing %s...'% strain;

  data[strain] = exp.load_data(strain);
  dat[strain] = getattr(data[strain], feat);
  rids = getattr(data[strain], 'roam') > 0;
  
  bins, bine = exp.stage_bins(data[strain], nbins= sbins);
  
  nworms = rids.shape[0];
  dat_bin[strain] = np.zeros((nworms, tbins));
  for t in range(tbins):
    for w in range(nworms):
      ids = np.arange(bins[w,t], bine[w,t]);
      dd = dat[strain][w, ids][rids[w][ids]];
      dat_bin[strain][w,t] = np.std(dd) / np.mean(dd);
    

  #dat_bin[strain] = exp.bin_data(dat[strain], sbinsb);
  dat_mean[strain] = np.mean(dat_bin[strain], axis = 0);
  dat_var[strain]  = np.var(dat_bin[strain], axis = 0);
  
  #stage_bins[strain] = exp.stage_bins(data[strain], nbins=1);
  #dat_bin_s[strain] = exp.bin_data(dat[strain], stage_bins[strain]);


#%% Order by activity

order = {};
for strain in strains:
  order[strain] = np.argsort(np.sum(dat[strain], axis = 1));
  
  
#%% Plot Data

# plot roaming speed
plt.figure(70); plt.clf()
fplt.plot_image_array([dat_bin[s][order[s]] for s in strains], names = strains, vmax = 10.0, cmap = 'jet')



#%%

for s in strains:
  print dat_bin[s].max()


#%% Calculate distributions

max_feat = 8.0;

#distributions
dist_t = {}; dist_s = {}; dist_g = {};
dist_a = np.zeros((dbins,1));

for strain in strains:
  #time resolved distribution
  dist_t[strain] = np.zeros((dbins, tbins));
  for t in range(tbins):
    dist_t[strain][:,t] = np.histogram(dat_bin[strain][:,t], range = (0, max_feat), bins = dbins)[0];
    dist_t[strain][:,t] /= np.sum(dist_t[strain][:,t]);
  
  #stage resolved distributions
  dist_s[strain] = np.zeros((dbins, nstages));
  for s in range(nstages):
    #  dist_s[strain][:,s] = np.histogram(dat_bin_s[strain][:,s], range = (0, max_feat), bins = dbins)[0];
    dist_s[strain][:,s] = np.histogram(dat_bin[strain][:,s*sbins:(s+1)*sbins].flatten(), range = (0, max_feat), bins = dbins)[0];
    dist_s[strain][:,s] /= np.sum(dist_s[strain][:,s]); 
  
  #genotype resolved distributions
  dist_g[strain] =  np.sum(dist_t[strain], axis = 1);
  dist_g[strain] =  np.array([dist_g[strain] / np.sum(dist_g[strain])]).T;
  
  dist_a += dist_g[strain];
dist_a /= np.sum(dist_a);
  


#%% Plot data on top of distributions


fig = plt.figure(6); plt.clf();

sts = np.array(strains)[1:]

for i, s in enumerate(sts):
  ax = plt.subplot(len(sts),1,i+1);
  im = plt.imshow(dist_t[s], vmin = 0, vmax = 1, origin = 'lower', aspect = 'auto', cmap = plt.cm.Oranges)
  # plot data on top
  xdat = np.array([np.ones(dat_bin[s].shape[0]) * t for t in range(dat_bin[s].shape[1])]).flatten()
  plt.scatter(xdat, dat_bin[s].T.flatten() * (dbins-0.5) /max_feat, c = 'gray', s = 2)
  plt.xlim(-0.5, tbins-0.5)
  plt.ylim(-0.5, dbins-0.5)
  
  plt.yticks(np.linspace(0,dbins,6)-0.5, np.linspace(0,max_feat,6));
  
  #if i < len(sts)-1:
  #plt.setp(ax.get_xticklabels(), visible=False)
  #plt.xticks(np.arange(0,tbins+1,)-0.5, ('', '','', '','', '', '', '', '', '', ''))
  #else:
  #plt.xticks(np.arange(0,51,5)-0.5, ('', 'L1','', 'L2','', 'L3', '', 'L4', '', 'A', ''))

  # add average trace
  plt.plot(np.mean(dat_bin[s] /max_feat *(dbins-0.5), axis = 0), 'k')
  #plt.title(s)
  plt.ylabel(s)
 
#plt.tight_layout()


# plot distributions 

#plt.figure(7); plt.clf();
#plt.plot(dist_t[strain][:,:10:20])

fig.subplots_adjust(right=0.5)
cbar_ax = fig.add_axes([0.09, 0.135, 0.22, 0.01])
fig.colorbar(im, cax=cbar_ax, orientation = 'horizontal')

fig.show()

#%%
fig.savefig(os.path.join(fig_dir, '%s_cv_roaming_speed_distributions_sbins=%d.pdf' % (strains_label, sbins)))


#%% Distance in time accross geno types to wild type 

nst = 0;
strain_ref = strains[nst]

cols = ['gray', 'k', 'lightblue', 'r', 'g', 'darkblue', 'black'];

strains_order = strains[1:];

#for nst, strain_ref in enumerate(strains):
for nst, strain_ref in enumerate([strains[0]]):
  ent_t_ref = {};
  N = {};
  for strain in strains:
    ent_t_ref[strain] = np.zeros(tbins);
    for t in range(tbins):
      ent_t_ref[strain][t] = d(dist_t[strain][:,t], dist_t[strain_ref][:,t]);
    N[strain] = dat_bin[strain].shape[0] + dat_bin[strain_ref].shape[0]
  
      
  # significance level form theory
      
  from scipy.special import gamma, gammainc, gammaincc
  
  # get counts for each time pairs
  sJS = {}; # the probability of observing this or a smaller value of dJS assuming distributions are not different
  nu = (2-1)*(dbins-1);
  for strain in strains:
    sJS[strain] = np.zeros(tbins);
    for t in range(tbins):
      sJS[strain][t] = gammainc(nu / 2.0, N[strain] * np.log(2) * np.abs(ent_t_ref[strain][t])); # / gamma(nu/2.0);
  
  
  p = 0.05;
  
  fig = plt.figure(200+10*nst); plt.clf();
  #plt.subplot(2,1,1)
  for ns, strain in enumerate(strains):
    plt.plot(ent_t_ref[strain], label = strain, c = cols[ns]);
    ids = np.where(1 - sJS[strain] <= p)[0];
    if len(ids) > 0:
      plt.plot(ids, ent_t_ref[strain][ids], 'o' ,c = cols[ns], markersize = 6)
    
  #plt.fill_between(range(tbins), 0, s_boot, color = 'lightgray')    
    
  plt.legend();
  plt.xlim(0, tbins);
  
  fig.savefig(os.path.join(fig_dir, '%s_cv_roaming_speed_JS_time_signed_ref=%s_sbins=%d.pdf' % (strains_label, strain_ref, sbins)));
  
  #plt.subplot(2,1,2);
  #plt.title('dJS - ref =%s' % strain_ref)
  #for ns, strain in enumerate(strains):
  #  plt.plot(1-sJS[strain], label = strain, c = cols[ns]);
  #plt.legend();
  #plt.title('p value - ref=%s' % strain_ref);
  
  
  # matrix plot
  fig = plt.figure(300); plt.clf();
  
  #ordr = np.array([5, 1, 2, 3, 4]); 
  #strains_order = [strains[i] for i in ordr];  
  dd = np.zeros((tbins, nstrains-1));
  for ns, strain in enumerate(strains_order):
    dd[:,ns] = ent_t_ref[strain];
    ids = np.where(1 - sJS[strain] > p)[0];
    if len(ids) > 0:
      dd[ids, ns] = 0;
  
  
  import matplotlib as mpl;
  
  wl = 0.85; wl2 = 0.5;
  f1 = 0.5; fl2 = 1.0;
  cdict1 = {'red':  ((0.0, 0.0, 0.0),
                     (0.25, wl2, wl2),
                     (0.5,  wl, wl),
                     (0.75, fl2, fl2),
                     (1.0,  f1, f1)),

            'green':((0.0, 0.0, 0.0),
                     (0.25, wl2, wl2),
                     (0.5,  wl, wl),
                     (0.75, wl2, wl2),
                     (1.0,  0.0, 0.0)),

           'blue':  ((0.0, f1, f1),
                     (0.25, fl2, fl2),
                     (0.5,  wl, wl),
                     (0.75, wl2, wl2),
                     (1.0,  0.0, 0.0))
        }
  
  clmp  = mpl.colors.LinearSegmentedColormap('cm', cdict1, N=256, gamma=1.0)  
  #clmp = plt.cm.coolwarm;
  
  plt.imshow(dd.T, aspect = 'auto', cmap = clmp, vmax = 1.0, vmin = -1.0)
  plt.tight_layout();
  plt.yticks(range(nstrains-1), strains_order)
  tks = np.linspace(-0.5,tbins-0.5, 6);
  for t in tks[1:-1]:
    plt.plot([t,t],[-0.5, nstrains-1.5], 'k')
  
  tks = (tks[1:] + tks[:-1])/2;
  plt.xticks(tks, ['L1', 'L2', 'L3', 'L4', 'Adult'])
  plt.colorbar()
  plt.xlim(-0.5, tbins-0.5);
  
  fig.savefig(os.path.join(fig_dir, '%s_cv_roaming_speed_JS_time_signed_ref=%s_matrix_sbins=%d.pdf' % (strains_label, strain_ref, sbins)));
  
  np.savetxt(os.path.join(fig_dir, '%s_cv_roaming_speed_JS_time_signed_ref=%s_matrix_sbins=%d.txt' % (strains_label, strain_ref, sbins)), dd)




  # matrix plot - significance
  fig = plt.figure(303); plt.clf();
  
  #ordr = np.array([5, 1, 2, 3, 4]); 
  #strains_order = [strains[i] for i in ordr];  
  dd = np.zeros((tbins, nstrains-1));
  for ns, strain in enumerate(strains_order):
    dd[:,ns] = (sJS[strain]) * np.sign(ent_t_ref[strain]);
    #ids = np.where(1 - sJS[strain] > p)[0];
    #if len(ids) > 0:
    #  dd[ids, ns] = 0;
  
  plt.imshow(dd.T, aspect = 'auto', cmap = plt.cm.coolwarm, vmax = 1.0, vmin = -1.0)
  plt.tight_layout();
  plt.yticks(range(nstrains-1), strains_order)
  tks = np.linspace(0,tbins, 6)- 0.5;
  for t in tks[1:-1]:
    plt.plot([t,t],[-0.5, nstrains-1.5], 'k')
  
  tks = (tks[1:] + tks[:-1])/2;
  plt.xticks(tks, ['L1', 'L2', 'L3', 'L4', 'Adult'])
  plt.colorbar()
  plt.xlim(-0.5, tbins-0.5);
  
  fig.savefig(os.path.join(fig_dir, '%s_cv_roaming_speed_JS_time_signed_significance_ref=%s_matrix_sbins=%d.pdf' % (strains_label, strain_ref, sbins)));
  


#%% Write data 

dd = np.zeros((len(strains_order),  ent_t_ref[strains_order[0]].shape[0]));
for i, s in enumerate(strains_order):
  dd[i] = ent_t_ref[s];
  
#from scipy.special import gammaincc;
pvals = np.zeros((len(strains_order), sJS[strains_order[0]].shape[0]));
for i,s in enumerate(strains_order):
  pvals[i] = gammaincc(nu / 2.0, N[s] * np.log(2) * np.abs(ent_t_ref[s]));
  
qvals = pvals.copy()
for i in range(qvals.shape[0]):
  qvals[i,:] = als.correctPValues(qvals[i]);
  
np.savetxt(os.path.join(fig_dir, '%s_cv_roaming_speed_JS_sbins=%d.txt' % (strains_label, sbins)), dd)
np.savetxt(os.path.join(fig_dir, '%s_cv_roaming_speed_JS_pvals_sbins=%d.txt' % (strains_label, sbins)), pvals)
np.savetxt(os.path.join(fig_dir, '%s_cv_roaming_speed_JS_qvals_sbins=%d.txt' % (strains_label, sbins)), qvals)
 
 

#%% Bottstrap

for nst, strain_ref in enumerate([strains[0]]):
  # significance via boot strapping
  # note: could start from roaming 0,1 s (but ok)
  
  nboots = 10000;
  dJS_boot = np.zeros((nboots, tbins));
  
  for t in range(tbins):
    print 'processing %d / %d' % (t, tbins);
    for b in range(nboots):
      dat_boot = np.random.choice(data_all[strain_ref][:,t], size = data_all[strain].shape[0], replace = True);
      dist_boot = np.array(np.histogram(dat_boot, range = (0, max_feat), bins = dbins)[0], dtype = float);
      dist_boot /= np.sum(dist_boot);
      dJS_boot[b,t] = np.abs(d(dist_boot, np.abs(dist_t[strain_ref][:,t])));
    
  plt.figure(201+10*nst); plt.clf();
  plt.plot(dJS_boot, '.')
  plt.title('bootstrapping - ref=%s' % strain_ref);
  
  # significance level:
  
  s_boot = np.percentile(dJS_boot, 90, axis = 0);
  
  plt.figure(202+10*nst); plt.clf();
  plt.fill_between(range(tbins), 0, s_boot, color = 'lightgray')
  for strain in strains:
    plt.plot(np.abs(ent_t_ref[strain]), label = strain, linewidth = 2);
    ids = np.where(np.abs(ent_t_ref[strain]) >= s_boot)[0];
    if len(ids) > 0:
      plt.plot(ids, np.abs(ent_t_ref[strain][ids]), '*k', markersize = 10)
  plt.legend();
  plt.xlim(-0.5, tbins-0.5)
  plt.title('dJS bootstrapping - ref=%s' % strain_ref);


