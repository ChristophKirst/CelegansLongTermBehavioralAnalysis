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


### Parameter

strains = exp.strains;
feat = 'roam';

nstrains = len(strains);
nstages = 5;

dbins = 20;
sbins = 10;
tbins = nstages * sbins;

fig_dir = '/data/Science/Projects/CElegansBehaviour/Analysis/DwellingRoaming/Figures/2016_12_20'
fig_dir = None;

#d = als.jensen_shannon_divergence;

def d(p,q):
  dd = als.jensen_shannon_divergence(p,q);
  if np.sum(np.linspace(0,1, len(p)) * p) > np.sum(np.linspace(0,1, len(q)) * q):
    return dd;
  else:
    return -dd;



### Load / Prepare data

data = {}; dat_bin = {}; dat_mean = {}; dat_var = {};
dat = {}; stage_bins = {}; dat_bin_s = {};

for strain in strains:
  print 'processing %s...'% strain;

  data[strain] = exp.load_data(strain);
  dat[strain] = getattr(data[strain], feat);
  
  sbinsb = exp.stage_bins(data[strain], nbins= sbins);

  dat_bin[strain] = exp.bin_data(dat[strain], sbinsb);
  dat_mean[strain] = np.mean(dat_bin[strain], axis = 0);
  dat_var[strain]  = np.var(dat_bin[strain], axis = 0);
  
  stage_bins[strain] = exp.stage_bins(data[strain], nbins=1);
  dat_bin_s[strain] = exp.bin_data(dat[strain], stage_bins[strain]);


### Order by activity

order = {};
for strain in strains:
  order[strain] = np.argsort(np.sum(dat[strain], axis = 1));
  
  
### Plot Data

# plot roaming fraction
plt.figure(70); plt.clf()
fplt.plot_image_array([dat_bin[s][order[s]] for s in strains], names = strains, vmax = 1.0, cmap = 'jet')




### Calculate distributions

max_feat = 1.0;

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
  

### Plot data on top of distributions

strain = strains[0];
#strain = strains[];

plt.figure(6); plt.clf();
sts = [strains[0],strains[-1]];
for i, s in enumerate(sts):
  plt.subplot(len(sts),1,i+1);
  plt.imshow(dist_t[s], vmax = 0.5, origin = 'lower', aspect = 'auto')
  # plot data on top

  plt.plot(dat_bin[s].T * (dbins-0.5), '.', c = 'w')
  plt.xlim(-0.5, tbins-0.5)
  plt.ylim(-0.5, dbins-0.5)

  # add average trace
  plt.plot(np.mean(dat_bin[s]*(dbins-0.5), axis = 0), 'r')
  plt.title(s)
  plt.tight_layout()



# plot distributions 

plt.figure(7); plt.clf();
plt.plot(dist_t[strain][:,:10:20])



### Distance in time accross geno types to wild type 


strain_ref = strains[0]

for nst, strain_ref in enumerate(strains):
  ent_t_ref = {};
  N = {};
  for strain in strains:
    ent_t_ref[strain] = np.zeros(tbins);
    for t in range(tbins):
      ent_t_ref[strain][t] = d(dist_t[strain][:,t], dist_t[strain_ref][:,t]);
    N[strain] = dat_bin[strain].shape[0] + dat_bin[strain_ref].shape[0]
  
      
  # significance level form theory
      
  from scipy.special import gamma, gammainc
  
  # get counts for each time pairs
  sJS = {}; # the probability of observing this or a smaller value of dJS assuming distributions are not different
  nu = (2-1)*(dbins-1);
  for strain in strains:
    sJS[strain] = np.zeros(tbins);
    for t in range(tbins):
      sJS[strain][t] = gammainc(nu / 2.0, N[strain] * np.log(2) * ent_t_ref[strain][t]); # / gamma(nu/2.0);
  
  
  p = 0.05;
  
  plt.figure(200+10*nst); plt.clf();
  plt.subplot(2,1,1)
  for strain in strains:
    plt.plot(ent_t_ref[strain], label = strain);
    ids = np.where(1 - sJS[strain] <= p)[0];
    if len(ids) > 0:
      plt.plot(ids, ent_t_ref[strain][ids], '*k', markersize = 10)
  plt.legend();
  plt.subplot(2,1,2);
  plt.title('dJS - ref =%s' % strain_ref)
  for strain in strains:
    plt.plot(1-sJS[strain], label = strain);
  plt.legend();
  plt.title('p value - ref=%s' % strain_ref);
  
  
  # significance via boot strapping
  # note: could start from roaming 0,1 s (but ok)
  
  nboots = 10000;
  dJS_boot = np.zeros((nboots, tbins));
  
  for t in range(tbins):
    print 'processing %d / %d' % (t, tbins);
    for b in range(nboots):
      dat_boot = np.random.choice(dat_bin[strain_ref][:,t], size = dat_bin[strain].shape[0], replace = True);
      dist_boot = np.array(np.histogram(dat_boot, range = (0, max_feat), bins = dbins)[0], dtype = float);
      dist_boot /= np.sum(dist_boot);
      dJS_boot[b,t] = d(dist_boot, dist_t[strain_ref][:,t]);
    
  plt.figure(201+10*nst); plt.clf();
  plt.plot(dJS_boot, '.')
  plt.title('bootstrapping - ref=%s' % strain_ref);
  
  # significance level:
  
  s_boot = np.percentile(dJS_boot, 99, axis = 0);
  
  plt.figure(202+10*nst); plt.clf();
  plt.fill_between(range(tbins), 0, s_boot, color = 'lightgray')
  for strain in strains:
    plt.plot(ent_t_ref[strain], label = strain, linewidth = 2);
    ids = np.where(ent_t_ref[strain] >= s_boot)[0];
    if len(ids) > 0:
      plt.plot(ids, ent_t_ref[strain][ids], '*k', markersize = 10)
  plt.legend();
  plt.xlim(-0.5, tbins-0.5)
  plt.title('dJS bootstrapping - ref=%s' % strain_ref);







# fit a beta distribution to distributions
import scipy.stats as stats

pdf_par = {};
for strain in strains:
  pdf_par[strain] = np.zeros((2, tbins));

  for t in range(tbins):
    fitdat = dat_bin[strain][:,t] + 10**-5;
    fitdat[fitdat > 1.0] = 1.0 - 10**-5;
    betafit = stats.beta.fit(fitdat,floc=0, fscale=1);
    pdf_par[strain][:,t] = [betafit[0], betafit[1]];


strain = strains[0]
plt.figure(8); plt.clf();
for t in range(tbins):
  plt.subplot(7,8,t+1)
  x = np.linspace(0,1,150);
  plt.plot(x, stats.beta.pdf(x, pdf_par[strain][0,t], pdf_par[strain][1,t], loc = 0, scale = 1))
  plt.plot(np.linspace(0,1,dbins), dist_t[strain][:,t] * dbins)


plt.figure(10); p:lt.clf();
for i,strain in enumerate(strains):
  plt.subplot(nstrains,2,1+i*2)
  plt.plot(pdf_par[strain][0,:])
  plt.plot(np.ones(tbins))
  plt.title('%s - alpha' % strain)
  plt.subplot(nstrains,2,2+i*2)
  plt.plot(pdf_par[strain][1,:])  
  plt.plot(np.ones(tbins))
  plt.title('%s - beta' % strain)


# pdf trjectories
plt.figure(11); plt.clf()
for strain in strains:
  plt.plot(pdf_par[strain][0,:],pdf_par[strain][1,:], label = strain)
plt.legend()
plt.xlabel('alpha')
plt.ylabel('beta')
  

#
plt.figure(11); plt.clf();
for i,strain in enumerate(strains):
  plt.subplot(1,2,1)
  plt.plot(pdf_par[strain][0,:], label = strain)
  plt.subplot(1,2,2)
  plt.plot(pdf_par[strain][1,:], label = strain)  


plt.subplot(1,2,1)
plt.legend()
plt.plot(np.ones(tbins))
plt.title('alpha')
plt.subplot(1,2,2)
plt.plot(np.ones(tbins))
plt.legend()
plt.title('alpha')



### distances to average roaming fraction

delta = {};
for strain in strains:
  delta[strain] = dat_bin[strain]-dat_mean[strain];

strain = strains[0];
plt.figure(60); plt.clf()
#plt.plot(delta[strain])
plt.hist(delta[strain].flatten(), bins = 64)

# z score in each time bin 
zscore = {}; mean = {}; var = {};
for strain in strains:
  mean[strain] =  np.mean(delta[strain], axis = 0);
  var[strain] = np.var(delta[strain], axis = 0);
  zscore[strain] = (delta[strain] - mean[strain])/ var[strain]


strain = strains[0];
plt.figure(111); plt.clf();
plt.subplot(1,3,1)
plt.plot(zscore[strain].T)
plt.ylim(-100,300)
plt.subplot(1,3,2)
plt.plot(mean[strain])
plt.subplot(1,3,3)
plt.plot(var[strain])


# clean lethragus etc
strain = strains[0];
var0 = var[strain][0]+10**-10;
ids = var[strain] <= var0;

zscore0 = {}; mean0 = {}; var0 = {};
for strain in strains:
  zscore0[strain] = zscore[strain].copy();
  zscore0[strain][:,ids] = 0;
  mean0[strain] = mean[strain].copy();
  mean0[strain][ids] = 0;
  var0[strain] = var[strain].copy();
  var0[strain][ids] = 0;
  

strain = strains[0];
plt.figure(112); plt.clf();
plt.subplot(1,3,1)
plt.plot(zscore0[strain].T)
plt.subplot(1,3,2)
plt.plot(mean0[strain])
plt.subplot(1,3,3)
plt.plot(var0[strain])


# see if we can find some consistency in z scores:
cons = {};
for strain in strains:
  cons[strain] = np.sum(zscore0[strain], axis = 1);

strain = strains[0]
plt.figure(113); plt.clf();

plt.subplot(1,2,1)
plt.plot(cons[strain], '.r')
plt.plot(np.zeros(len(cons[strain])), 'k')
plt.xlim(0, len(cons[strain]))

plt.subplot(1,2,2); 
plt.hist(cons[strain])


# boostrapped consistency index !

strain = strains[0];

k = 0;
for strain in strains:
  
  n_boot_strap = 50000;
  sc = zscore0[strain].copy();
  sc = sc[:,~ids];
  sc = sc.flatten();
  ns = zscore0[strain].shape[1] - ids.sum();
  
  bs = np.random.choice(sc, size = (n_boot_strap, ns), replace = True);
  
  cons_bs = {};
  cons_bs[strain] = np.sum(bs, axis = 1);
  
  nw = len(cons[strain]);
  

  
  hbins = 10;
  hrange = [-300, 500];
  
  plt.figure(113+k); plt.clf();
  
  plt.subplot(1,2,1)
  plt.plot(np.linspace(0, nw, n_boot_strap), cons_bs[strain], '.', c='gray')
  plt.plot(np.zeros(nw), 'k')
  plt.xlim(0, len(cons[strain]))
  
  plt.subplot(1,2,2); 
  plt.hist(cons_bs[strain], bins = hbins, normed = True, color = 'gray', alpha = 0.5, range = hrange)
  
  
  plt.subplot(1,2,1)
  plt.plot(cons[strain], '.r', markersize= 20)
  plt.plot(np.zeros(len(cons[strain])), 'k')
  plt.xlim(0, len(cons[strain]))
  
  plt.subplot(1,2,2); 
  plt.hist(cons[strain], bins = hbins, normed = True, color = 'red', alpha = 0.5,  range = hrange)
  
  
  # calculate the significant outliers
  
  pc = np.percentile(cons_bs[strain], [1, 5, 95, 99])
  
  plt.subplot(1,2,1);
  for p in pc:
    plt.plot(p * np.ones(nw), 'b');
  
  
  # difference of distribution for worms to bootstrapped one ?
  
  hw = np.histogram(cons[strain], bins = hbins, normed = True, range = hrange)[0]
  hb = np.histogram(cons_bs[strain], bins = hbins, normed = True, range = hrange)[0]
  
  dd = d(hw, hb)
  
  
  nu = (2-1)*(hbins-1);
  ssJS = gammainc(nu / 2.0, (nw * 2) * np.log(2) * dd);
  
  plt.title(strain)
  plt.subplot(1,2,2);
  plt.title('d JS: %f; significance level %f' % (dd, ssJS))
  k += 1;
  
  


#fplt.plot_trace(pdf_par[strain].T)

#distances  
ent_t = {}; ent_s = {}; ent_g = {}; 
ent_t_s = {}; ent_t_g = {}; ent_s_g = {}
ent_t_a = {}; ent_s_a = {}; ent_g_a = {}
ent_t_0 = {};


dist_0 = np.zeros(dbins);
dist_0[0] = 1.0;

for strain in strains:
  print 'processing %s...' % strain
  ent_t[strain] = np.zeros((tbins, tbins));
  for t in range(tbins):
    for t2 in range(tbins):
       ent_t[strain][t,t2] = d(dist_t[strain][:,t], dist_t[strain][:,t2]);

  ent_s[strain] = np.zeros((nstages, nstages));
  for s in range(nstages):
    for s2 in range(nstages):
       ent_s[strain][s,s2] = d(dist_s[strain][:,s], dist_s[strain][:,s2]);
  
  ent_t_s[strain] = np.zeros((tbins, nstages))
  for t in range(tbins):
    for s in range(nstages):
      ent_t_s[strain][t,s] = d(dist_t[strain][:,t], dist_s[strain][:,s]);
    
  ent_t_g[strain] = np.array([d(dist_t[strain][:,t], dist_g[strain][:,0]) for t in range(tbins)]);
  ent_s_g[strain] = np.array([d(dist_s[strain][:,s], dist_g[strain][:,0]) for s in range(nstages)]);
  
  ent_t_a[strain] = np.array([d(dist_t[strain][:,t], dist_a[:,0]) for t in range(tbins)]);
  ent_s_a[strain] = np.array([d(dist_s[strain][:,s], dist_a[:,0]) for s in range(nstages)]);
  
  ent_g_a[strain] = d(dist_g[strain], dist_a[:,0]);
  
  ent_t_0[strain] = np.array([d(dist_t[strain][:,t], dist_0) for t in range(tbins)]);

  


#distances to wild type
sref = 'N2';
ent_t_ref = {}; ent_s_ref = {}; ent_g_ref = {}; 
ent_t_s_ref = {}; ent_t_g_ref = {}; ent_s_g_ref = {}

for strain in strains:
  print 'processing %s...' % strain
  ent_t_ref[strain] = np.zeros((tbins, tbins));
  for t in range(tbins):
    for t2 in range(tbins):
       ent_t_ref[strain][t,t2] = d(dist_t[strain][:,t], dist_t[sref][:,t2]);

  ent_s_ref[strain] = np.zeros((nstages, nstages));
  for s in range(nstages):
    for s2 in range(nstages):
       ent_s_ref[strain][s,s2] = d(dist_s[strain][:,s], dist_s[sref][:,s2]);
  
  ent_t_s_ref[strain] = np.zeros((tbins, nstages))
  for t in range(tbins):
    for s in range(nstages):
      ent_t_s_ref[strain][t,s] = d(dist_t[strain][:,t], dist_s[sref][:,s]);
    
  ent_t_g_ref[strain] = np.array([d(dist_t[strain][:,t], dist_g[sref][:,0]) for t in range(tbins)]);
  ent_s_g_ref[strain] = np.array([d(dist_s[strain][:,s], dist_g[sref][:,0]) for s in range(nstages)]);
  
  

#full distance matrix
def full_distance(dist_f, fbins):
  fbins_a = len(strains) * fbins;
  ent = np.zeros((fbins_a,fbins_a));
  for i,s1 in enumerate(strains):
    for j,s2 in enumerate(strains):
      print '%s %s' % (s1,s2);
      for t in range(fbins):
        for t2 in range(fbins):
          if j < i:
            ent[t+i*fbins,t2 + j*fbins] = d(dist_f[s1][:,t], dist_f[s2][:,t2]);
          elif i == j:
            if t2 <= t:
              ent[t+i*fbins,t2 + j*fbins] = d(dist_f[s1][:,t], dist_f[s2][:,t2]);

  irows,icols = np.triu_indices(len(ent),1)
  ent[irows,icols]=ent[icols,irows];
  return ent;

ent_t['all'] = full_distance(dist_t, tbins);
ent_s['all'] = full_distance(dist_s, nstages);
ent_g['all'] = full_distance(dist_g, 1);



# Plot densities and entropy differences


fig = plt.figure(70); plt.clf();
for strain in strains:
  plt.plot(ent_t_0[strain])







# Plot Clustering at different time resolutions
fig = plt.figure(1); plt.clf();
ents = [ent_t, ent_s, ent_g];
for i,e in enumerate(ents):
  plt.subplot(1,len(ents),i+1);
  plt.imshow(1-np.sqrt(e['all']), interpolation = 'none', aspect = 'auto', cmap = plt.cm.viridis)
  plt.title('%r' % strains)
plt.tight_layout()




reload(fplt)
scales = ['time', 'stage', 'life']
for ei,e in enumerate(ents):
  mat = e['all'];
  sub_size = len(mat)/ nstrains;
  mat = np.split(mat, nstrains);
  mat = [np.split(m, nstrains, axis = 1) for m in mat];
  dd = np.array([[np.linalg.norm(np.sqrt(mat[i][j])) for i in range(nstrains)] for j in range(nstrains)]);
  dd = dd[np.triu_indices(len(dd), 1)];
  
  fig = plt.figure(2+ei); plt.clf();
  fplt.plot_hierarchical_cluster(e['all'], clustering = dd,  
                                 linkage_kwargs = {'method' : 'average', 'metric' : 'correlation'},
                                 label = exp.strains, colorbar_width = 0.02, padding = [0.05, 0.05], stride_line_kwargs={'c' : 'k', 'linewidth' : 2})
                                 
  fig.suptitle('JS diverence, clustering = || JSD[i,j] ||, scale = %s' % scales[ei])
  if fig_dir is not None:
    fig.savefig(os.path.join(fig_dir, 'JSD_norm_clustering_tbins=%d_scale=%s.pdf' % (tbins, scales[ei])))


reload(fplt)
for ei,e in enumerate(ents):
  mat = e['all'];
  sub_size = len(mat)/ nstrains;
  mat = np.split(mat, nstrains);
  mat = [np.split(m, nstrains, axis = 1) for m in mat];
  dd = np.sqrt(np.vstack([mat[i][i].flatten() for i in range(nstrains)]))

  #dd = dd[np.triu_indices(len(dd))];
  
  fig = plt.figure(5+ei); plt.clf();
  fplt.plot_hierarchical_cluster(e['all'], clustering = dd,  
                                 linkage_kwargs = {'method' : 'single', 'metric' : 'correlation'},
                                 label = exp.strains, colorbar_width = 0.02, padding = [0.05, 0.05], stride_line_kwargs={'c' : 'k', 'linewidth' : 2})
  
  fig.suptitle('JS diverence, clustering = JSD correlation, scale = %s' % scales[ei])
  if fig_dir is not None:
    fig.savefig(os.path.join(fig_dir, 'JSD_correlation_clustering_tbins=%d_scale=%s.pdf' % (tbins, scales[ei])))


### Cluster all distribution via similarity

fig = plt.figure(80); plt.clf();
dd = ent_t['all'].copy();
dd = dd[np.triu_indices(len(dd), 1)];

fplt.plot_hierarchical_cluster(ent_t['all'],clustering = dd,  
                                 linkage_kwargs = {'method' : 'average', 'metric' : 'correlation'},
                                 colorbar_width = 0.02, padding = [0.05, 0.05], stride_line_kwargs={'c' : 'k', 'linewidth' : 2})

fig.suptitle('JS diverence, clustering = all times')
if fig_dir is not None:
  fig.savefig(os.path.join(fig_dir, 'JSD_time_clustering_tbins=%d.pdf' % (tbins)))


### Manifold embeddings

reload(fplt);

plt.figure(200); plt.clf();
#label = range(tbins) * nstrains
#label = list(np.sort(range(nstages) * sbins)) * nstrains
label = np.sort(range(nstrains) * tbins)
res = fplt.plot_manifold_embeddings(ent_t['all'], precomputed=True, label = label)


reload(fplt);
fig = plt.figure(50); plt.clf();
label = np.sort(range(nstrains) * tbins)
colors = ['red', 'blue', 'darkgreen', 'orange', 'purple', 'gray'];
cmaps = [plt.cm.Reds, plt.cm.Blues, plt.cm.Greens, plt.cm.Oranges, plt.cm.Purples, plt.cm.gray_r];
i = 0;
for k,v, in res.iteritems():
  plt.subplot(1,3,i+1); i+=1;
  xmin = np.min(v[:,0]); xmax = np.max(v[:,0]);
  ymin = np.min(v[:,1]); ymax = np.max(v[:,1]);
  for l in range(nstrains):
    idx = label == l;
    fplt.plot_embedding(v[idx], title = k, color = colors[l]);
    fplt.plot_embedding_contours(v[idx], cmap = cmaps[l], contours = 15,  xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax);
#fig.suptitle('Embeddings JS-Diverence all')
if fig_dir is not None:
  fig.savefig(os.path.join(fig_dir, 'new_JSD_embedding_all_tbins=%d.pdf' % (tbins)))

  fig.savefig(os.path.join(fig_dir, 'new_JSD_embedding_all_tbins=%d.png' % (tbins)))

##full distributions
reload(fplt);
fig = plt.figure(50); plt.clf();
label = np.sort(range(nstrains) * tbins)
colors = ['red', 'blue', 'darkgreen', 'orange', 'purple', 'gray'];
cmaps = [plt.cm.Reds, plt.cm.Blues, plt.cm.Greens, plt.cm.Oranges, plt.cm.Purples, plt.cm.gray_r];
i = 0;
for k,v, in res.iteritems():
  plt.subplot(1,3,i+1); i+=1;
  xmin = np.min(v[:,0]); xmax = np.max(v[:,0]);
  ymin = np.min(v[:,1]); ymax = np.max(v[:,1]);
  #idx = label == l;
  fplt.plot_embedding(v, title = k, color = colors[l]);
  fplt.plot_embedding_contours(v, cmap = 'gray', contours = 15,  xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax);
#fig.suptitle('Embeddings JS-Diverence all')
if fig_dir is not None:
  fig.savefig(os.path.join(fig_dir, 'new_JSD_embedding_all_tbins=%d_combined.png' % (tbins)))



reload(fplt)
fig = plt.figure(51); plt.clf();
i = 0;
for k,v, in res.iteritems():
  xmin = np.min(v[:,0]); xmax = np.max(v[:,0]);
  ymin = np.min(v[:,1]); ymax = np.max(v[:,1]);
  for l in range(nstrains):
    plt.subplot(3,6,i+1); i+=1;
    idx = label == l;
    fplt.plot_embedding_contours(v[idx], cmap = cmaps[l], contours = None, density = True, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax);
    fplt.plot_colored_line(v[idx,0], v[idx,1], range(tbins), cmap = plt.cm.gray)
    fplt.plot_embedding(v[idx], color = 'gray');
    plt.title('%s %s' % (strains[l], k))
if fig_dir is not None:
  fig.savefig(os.path.join(fig_dir, 'new_JSD_embedding_all_strains_tbins=%d.png' % (tbins)))


reload(fplt)
fig = plt.figure(51); plt.clf();
i = 0;
for k,v, in res.iteritems():
  xmin = np.min(v[:,0]); xmax = np.max(v[:,0]);
  ymin = np.min(v[:,1]); ymax = np.max(v[:,1]);
  for l in range(nstrains):
    plt.subplot(3,6,i+1); i+=1;
    idx = label == l;
    fplt.plot_embedding_contours(v[idx], cmap = cmaps[l], contours = None, density = True);
    fplt.plot_colored_line(v[idx,0], v[idx,1], range(tbins), cmap = plt.cm.gray)
    fplt.plot_embedding(v[idx], color = 'gray');
    plt.title('%s %s' % (strains[l], k))
if fig_dir is not None:
  fig.savefig(os.path.join(fig_dir, 'new2_JSD_embedding_all_strains_tbins=%d.png' % (tbins)))


### Plot Distributions and distances

#for i, strain in enumerate(strains):

for i,s0 in enumerate(strains):
  fig = plt.figure(300+i, figsize = (21,14)); plt.clf();
  ax = plt.subplot(3,1,1);
  rs = plt.cm.Reds(0.9);
  fplt.plot_distributions(dat_bin[s0], percentiles_colors=['black', rs, 'red',  rs, 'black'],percentiles = [5,25,50,75,95], cmap = plt.cm.viridis)
  plt.xlim(0,tbins);
  plt.ylim(0,1);
  plt.title('%s %s' % (s0, feat))
  plt.subplot(3,1,2, sharex = ax)
  plt.plot(ent_t_g[s0])
  plt.title('JS divergence to average distribution')
  plt.xlim(0, len(ent_t_g[s0]))
  plt.subplot(3,1,3, sharex = ax)
  plt.plot(ent_t_s_ref[s0])
  plt.title('JS divergence to average stage distributions')
  plt.xlim(0, len(ent_t_g[s0]))
  plt.tight_layout();
    
  if fig_dir is not None:
    fig.savefig(os.path.join(fig_dir, 'Distribution_tbins=%d_strain=%s.pdf' % (tbins, s0)))

  fig = plt.figure(400+i, figsize = (21,14)); plt.clf();
  plt.plot(dist_s[s0])
  plt.plot(dist_g[s0], color = 'k', linewidth = 2)
  plt.title('Average distributions %s' % s0)
  plt.tight_layout()

  if fig_dir is not None:
    fig.savefig(os.path.join(fig_dir, 'Distribution_average_tbins=%d_strain=%s.pdf' % (tbins, s0)))






### Embedding + watershed classification
 
import scipy.stats as st

#print xmin,xmax,ymin,ymax

Y = res['t-SNE'];

npts = 500;
sf = 1.05;
xmin = sf*np.min(Y[:,0]); xmax = sf*np.max(Y[:,0]);
ymin = sf*np.min(Y[:,1]); ymax = sf*np.max(Y[:,1]);
dx = float(xmax-xmin) / npts;
dy = float(ymax-ymin) / npts;
xx, yy = np.mgrid[xmin:xmax:dx, ymin:ymax:dy]
positions = np.vstack([xx.ravel(), yy.ravel()])
#kernel = st.gaussian_kde(Y.T, bw_method = 'silverman');
kernel = st.gaussian_kde(Y.T, bw_method = 0.25);

Y_xy = kernel(positions)
Y_xy = np.reshape(Y_xy, xx.shape)

fplt.plot_array(Y_xy)

#watershed it
from scipy import ndimage as ndi
from skimage.morphology import watershed
from skimage.feature import peak_local_max


local_maxi = peak_local_max(Y_xy, indices=False, footprint=np.ones((7, 7)));
markers = ndi.label(local_maxi)[0]
labels = watershed(-Y_xy, markers) #, mask=image)

#classify points

Y_idx = np.array(np.round((Y - [xmin,ymin]) / [xmax-xmin, ymax-ymin] * (npts-1)), dtype = int);
Y_idx[Y_idx < 0] = 0; Y_idx[Y_idx > npts-1] = npts-1
Y_class = labels[Y_idx[:,0], Y_idx[:,1]]

cmap = 'jet'
fig = plt.figure(1); plt.clf();
plt.subplot(1,2,1);
plt.imshow(1-Y_xy.T, cmap = 'gray', extent=(xmin,xmax,ymin,ymax), origin = 'lower');
plt.xlim(xmin,xmax); plt.ylim(ymin,ymax)
plt.scatter(Y[:,0], Y[:,1], c= Y_class,  s = 30, edgecolors = 'face', cmap = cmap);
plt.subplot(1,2,2);
plt.imshow(labels.T,  extent=(xmin,xmax,ymin,ymax), origin = 'lower', cmap = cmap);


fig.savefig(os.path.join(fig_dir, 'cluster_2.png'))





## plot time evolution of the trajectories for each geno type

def Y_to_data(Y, nstrains, nbins):
  d = np.zeros((nstrains, nbins));
  for i in range(nstrains):
    for t in range(nbins):
      d[i,t] = Y[i*nbins + t];
  return d;

Y_dat = Y_to_data(Y_class, nstrains, tbins)


fig = plt.figure(5);
plt.imshow(Y_dat, cmap = cmap, aspect = 'auto')
plt.yticks(range(6), strains)



fig.savefig(os.path.join(fig_dir, 'cluster_2_in_time.png'))


































### Stretching of distirubtions


dist_sinh = {};    
    
    
  



### Calculate distances between distributions




### Plot Distributions and and distances















### Adaptive bin size distributions and Kernel density estimates 
# see: Freedman–Diaconis rule binsize = 2 IQR / n^(1/3)

iqr = {};
for strain in strains:
  iqr[strain] = np.percentile(dat_bin[strain], [5, 25, 50, 75, 95], axis = 0);


nworms = dat_bin[strains[0]].shape[0];
bin_width = {};
bin_width_max = {};
nbins_iqr = {};
nbins_iqr_max = {};
dist_iqr = {};
dist_iqr_max = {};
for strain in strains:
  bin_width[strain] = 2 * (iqr[strain][3,:]- iqr[strain][1,:]) / nworms**(1.0/3);
  bin_width_max[strain] = np.max(bin_width[strain]);
  bwm = bin_width[strain];
  bwm[bwm == 0] = np.max(dat_bin[strain]);
  nbins_iqr[strain] = np.ceil(np.max(dat_bin[strain], axis = 0) / bwm);
  nbins_iqr_max[strain] = np.max(nbins_iqr[strain]);

dbins_max = max([nbins_iqr_max[strain] for strain in strains]);
  
  
s0 = 'N2';
plt.figure(9); plt.clf();
plt.subplot(2,1,1);
plt.plot(bin_width[s0])
plt.subplot(2,1,2);
plt.plot(nbins_iqr[s0])
  
  
for strain in strains:  
  dist_iqr[strain] = np.zeros((nbins_max[strain], nbins));
  dist_iqr_max[strain] = np.zeros((nbins_max[strain], nbins));
  
  for b in range(nbins):
    dist_iqr[strain][:,b] += np.histogram(dat_bin[strain][:,b], range = (0, 0.9), bins = dbins)[0];
    dist_iqr[strain][:,b] /= np.sum(dist_iqr[strain][:,b]);



plt.figure(90); plt.clf();
s0 = 'N2'
fplt.plot_distributions(dat_bin[s0])
plt.ylim(0,1);
plt.xlim(0, nbins)
plt.tight_layout()
plt.plot(dat_mean[s0], c = 'k', linewidth = 2)



plt.figure(91); plt.clf();
s0 = 'N2'
fplt.plot_distributions(dat_bin_sinh, cmap = plt.cm.plasma)
plt.ylim(0,1);
plt.xlim(0, nbins)
plt.tight_layout()
plt.plot(dat_mean[s0], c = 'k', linewidth = 3)
plt.colorbar(fraction = 0.01, pad = 0.01)



from scipy.stats import gaussian_kde

kde = gaussian_kde(dat_bin[s0][:,100])(np.linspace(0,1,100));
plt.figure(92); plt.clf();
plt.plot(np.linspace(0,1,100), kde);
plt.scatter(dat_bin[s0][:,100], np.zeros(nworms), s = 30, c = 'k')





# simply plot all the points
plt.figure(100); plt.clf()
plt.subplot(3,1,1);
s0 = 'N2'
nworms, ntimes = dat_bin[s0].shape;
for i in range(ntimes):
  plt.scatter(i * np.ones(nworms), np.sort(dat_bin[s0][:,i]), c = range(nworms), cmap = plt.cm.Spectral, edgecolor = 'face')
for i in range(5):
  #plt.plot(iqr[s0][i,:],  c = plt.cm.Spectral(i/5.0), linewidth = 3);
  cols = ['gray', 'gray', 'red', 'gray', 'gray'];
  plt.plot(iqr[s0][i,:],  c = cols[i], linewidth = 3);
plt.xlim(0,ntimes);
plt.ylim(0,1);
plt.tight_layout();

plt.subplot(3,1,2);
eps = 0.1;
dat_bin_sinh = np.arcsinh(dat_bin[s0] / eps) / np.arcsinh(1.0/ eps);
iqr_sinh = np.percentile(dat_bin_sinh, [5, 25, 50, 75, 95], axis = 0);
for i in range(ntimes):
  plt.scatter(i * np.ones(nworms), np.sort(dat_bin_sinh[:,i]), c = range(nworms), cmap = plt.cm.Spectral, edgecolor = 'face')
for i in range(5):
  #plt.plot(iqr[s0][i,:],  c = plt.cm.Spectral(i/5.0), linewidth = 3);
  cols = ['gray', 'gray', 'red', 'gray', 'gray'];
  plt.plot(iqr_sinh[i,:],  c = cols[i], linewidth = 3);
plt.xlim(0,ntimes);
plt.ylim(0,1);
plt.tight_layout();

plt.subplot(3,1,3);
plt.imshow(dist[s0], interpolation = 'none', cmap = 'viridis', aspect = 'auto', origin = 'lower');






## MDS on the distance matrix
import sklearn.manifold as sm;
dbins = 10
dist_sinh = np.zeros((dbins, nbins));
for b in range(nbins):
  dist_sinh[:,b] += np.histogram(dat_bin_sinh[:,b], range = (0, 1.0), bins = dbins)[0];
  dist_sinh[:,b] /= np.sum(dist_sinh[:,b]);



ent_diff_dict = {}
for i,s1 in enumerate(['N2']):
  ent_diff = np.zeros((nbins,nbins));
  for b in range(nbins):
    for b2 in range(nbins):
      #ent_diff[b,b2] = d(dist[s1][:,b], dist[s2][:,b2]);
      ent_diff[b,b2] = d(dist_sinh[:,b], dist_sinh[:,b2]);
  
  ent_diff_dict[s1] = ent_diff;
  
  n_components = 2;
  mds = sm.MDS(n_components = n_components, max_iter=100, n_init=1, dissimilarity = 'precomputed')
  
  Y = mds.fit_transform(ent_diff)
  color = range(len(Y[:,0]));
  cmap = plt.cm.Spectral;
  fig = plt.figure(130+i); plt.clf();
  if n_components == 2:
    plt.subplot(1,2,1);
    plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=cmap, s= 100)
    plt.plot(Y[:, 0], Y[:, 1], 'k')
    plt.colorbar(pad = 0.1, fraction = 0.1);
  else:
    ax = fig.add_subplot(1,2,1, projection = '3d');
    ax.scatter(Y[:, 0], Y[:, 1], Y[:,2], c = color, cmap=cmap)
    plt.plot(Y[:, 0], Y[:, 1], Y[:,2], 'k') 
  plt.title("MDS %s %s" % (s1,s2))
  
  plt.subplot(1,2,2);
  plt.imshow(ent_diff, interpolation = 'none', cmap = plt.cm.magma, aspect = 'auto')
  plt.tight_layout();






fplt.plot_manifold_embeddings(ent_diff, dissimilarity='precomputed')

## try irq based resolution of bins






       

# difference betwen strain distributions

pairs = [['N2', 'tph-1'], ['N2', 'cat-2']]


for i,s12 in enumerate(pairs):
  s1, s2 = s12;

  ent_diff = np.zeros((nbins));
  for b in range(nbins):
    ent_diff[b] = d(dist[s1][:,b], dist[s2][:,b]);
  
  plt.figure(i+1); plt.clf();
  plt.subplot(5,1,1)
  plt.imshow(dat_bin[s1], interpolation = 'none', aspect = 'auto', cmap = 'viridis');
  plt.title(s1)
  plt.colorbar(pad = 0.01, fraction = 0.01);
  plt.subplot(5,1,2)
  plt.imshow(dat_bin[s2], interpolation = 'none', aspect = 'auto', cmap = 'viridis');
  plt.title(s2)
  plt.colorbar(pad = 0.01, fraction = 0.01);
  plt.subplot(5,1,3)
  plt.imshow(dist[s1], interpolation = 'none', aspect = 'auto', cmap = 'viridis', origin = 'lower');
  plt.title(s1)
  plt.colorbar(pad = 0.01, fraction = 0.01);
  plt.subplot(5,1,4)
  plt.imshow(dist[s2], interpolation = 'none', aspect = 'auto', cmap = 'viridis', origin = 'lower');
  plt.title(s2)
  plt.colorbar(pad = 0.01, fraction = 0.01);
  plt.subplot(5,1,5)
  plt.plot(ent_diff);
  plt.title('JS divergence')
  plt.ylim(0, 0.6);
  plt.xlim(0, len(ent_diff)-1)
  plt.tight_layout();





## MDS on the distance matrix
import sklearn.manifold as sm;
  
strains = exp.strains;

ent_diff_dict = {}

for i,p in enumerate(strains):
  s1 = s2 = p;
  ent_diff = np.zeros((nbins,nbins));
  for b in range(nbins):
    for b2 in range(nbins):
      ent_diff[b,b2] = d(dist[s1][:,b], dist[s2][:,b2]);
  
  ent_diff_dict[tuple(p)] = ent_diff;
  
  n_components = 2;
  mds = sm.MDS(n_components, max_iter=100, n_init=1, dissimilarity = 'precomputed')
  
  Y = mds.fit_transform(ent_diff)
  color = range(len(Y[:,0]));
  cmap = plt.cm.Spectral;
  fig = plt.figure(130+i); plt.clf();
  if n_components == 2:
    plt.scatter(Y[:, 0], Y[:, 1], c=color, cmap=cmap, s= 100)
    plt.plot(Y[:, 0], Y[:, 1], 'k')
    plt.colorbar(pad = 0.1, fraction = 0.1);
  else:
    ax = fig.add_subplot(1,1,1, projection = '3d');
    ax.scatter(Y[:, 0], Y[:, 1], Y[:,2], c = color, cmap=cmap)
  plt.title("MDS %s %s" % (s1,s2))
  plt.tight_layout();



## MDS on the full distance matrix and corresponding trajectories in those spaces
import sklearn.manifold as sm;
  
full = exp.strains;

# calc distance metrics

nbins_full = len(full) * nbins;
ent_diff_full = np.zeros((nbins_full,nbins_full));
for i,s1 in enumerate(full):
  for j,s2 in enumerate(full):
    for b in range(nbins):
      for b2 in range(nbins):
        ent_diff_full[b+i*nbins,b2 + j*nbins] = d(dist[s1][:,b], dist[s2][:,b2]);
  
 
n_components = 2;
mds = sm.MDS(n_components, max_iter=100, n_init=1, dissimilarity = 'precomputed')

Y_full = mds.fit_transform(ent_diff_full)

color = range(nbins);
cmap = plt.cm.Spectral_r;
fig = plt.figure(200); plt.clf();

for i,s1 in enumerate(full):
  if n_components == 2:
    plt.subplot(3,2,i+1);
    plt.scatter(Y_full[(i*nbins):((i+1)*nbins), 0], Y_full[(i*nbins):((i+1)*nbins), 1], c=color, cmap=cmap, s= 100)
    plt.plot(Y_full[(i*nbins):((i+1)*nbins), 0], Y_full[(i*nbins):((i+1)*nbins), 1], 'k')
    plt.colorbar(pad = 0.1, fraction = 0.1);
  else:
    ax = fig.add_subplot(3,2,i+1, projection = '3d');
    ax.scatter(Y_full[(i*nbins):((i+1)*nbins), 0], Y_full[(i*nbins):((i+1)*nbins), 1], Y_full[(i*nbins):((i+1)*nbins),2], c = color, cmap=cmap)
  plt.title("MDS %s" % s1);
  plt.tight_layout();

# plot all points lable by stage

color = range(nbins);
cmap = plt.cm.Spectral_r;
fig = plt.figure(200); plt.clf();

slabel = np.repeat(range(5), nbins/7);
slabel_full = np.array(list(slabel) * 5)

for i,s1 in enumerate(full):
  if n_components == 2:
    plt.subplot(3,2,i+1);
    plt.scatter(Y_full[(i*nbins):((i+1)*nbins), 0], Y_full[(i*nbins):((i+1)*nbins), 1], c=color, cmap=cmap, s= 100)
    plt.plot(Y_full[(i*nbins):((i+1)*nbins), 0], Y_full[(i*nbins):((i+1)*nbins), 1], 'k')
    plt.colorbar(pad = 0.1, fraction = 0.1);
    plt.xlim(-0.3, 0.5);
    plt.ylim(-0.5, 0.3);
  else:
    ax = fig.add_subplot(3,2,i+1, projection = '3d');
    ax.scatter(Y_full[(i*nbins):((i+1)*nbins), 0], Y_full[(i*nbins):((i+1)*nbins), 1], Y_full[(i*nbins):((i+1)*nbins),2], c = color, cmap=cmap)
  plt.title("MDS %s" % s1);
  plt.tight_layout();




fplt.plot_pca()






plt.figure(6); plt.clf();
fplt.plot_image_array((dist,), names = ['dist_%s' % feat], invert_y = True, nplots = 3);
plt.subplot(3,1,2);
plt.plot(dat_mean, 'r');
plt.title('%s_%s mean' % (strain, feat))
plt.subplot(3,1,3)
plt.plot(dat_var, 'b');
plt.title('%s_%s var' % (strain, feat))


# entropies /cross entropies

ent = np.zeros((nbins, nbins));
for b in range(nbins):
  for b2 in range(nbins):
    if b==b2:
      ent[b,b2] = stats.entropy(dist[:,b]);
    else:
      ent[b,b2] = stats.entropy(dist[:,b], dist[:,b2]);


plt.figure(7); plt.clf();
fplt.plot_image_array((dist, ent), names = ['dist_roam', 'entropy'], invert_y = True, nplots = 5);
plt.subplot(5,1,3);
plt.plot(dat_mean, 'r');
plt.title('%s_%s mean' % (strain, feat))
plt.subplot(5,1,4)
plt.plot(dat_var, 'b');
plt.title('%s_%s var' % (strain, feat))
plt.subplot(5,1,5)
plt.plot(np.diag(ent), 'b');
plt.title('%s_%s ent' % (strain, feat))


  


plt.figure(17); plt.clf();
fplt.plot_image_array((dat_bin, dist), names = ['dat_bin', 'dist_roam'], invert_y = True, nplots = 6);
plt.subplot(6,1,3);
plt.plot(dat_mean, 'r');
plt.title('%s_%s mean' % (strain, feat))
plt.subplot(6,1,4)
plt.plot(np.diag(ent), 'b');
plt.title('%s_%s entropy' % (strain, feat))
plt.subplot(6,1,5)
plt.plot(np.append([0], np.diag(ent2, k=1)), 'b');
plt.title('%s_%s KL n=1' % (strain, feat))
plt.subplot(6,1,6)
plt.plot(np.append([0,0], np.diag(ent2, k = 2)), 'b');
plt.title('%s_%s KL n=2' % (strain, feat))
plt.subplot(6,1,2)
dent = (np.append([0], np.diag(ent2, k = 1))); dent /= dent.max();

dent = np.convolve(dent, np.ones(5)/5.0, mode = 'same')

plt.plot(dent * (dbins-1), 'w');
mm = dat_mean / 0.9 * (dbins-1); 
plt.plot(mm, 'r');
plt.xlim(0, len(dat_mean)-1);
plt.ylim(0, dbins-1);
plt.tight_layout();    
      

# plot all distribution / cluster distributions !!

plt.figure(8)
for b in range(nbins):
  plt.plot(dist[:,b])

C,T,pvar  = als.pca(dist)


plt.figure(9); plt.clf();
plt.subplot(2,2,1);
plt.imshow(C, interpolation = 'none', aspect = 'auto', cmap = 'viridis');
plt.title('spatial principal components')
plt.colorbar(pad = 0.01, fraction = 0.01);
ax = plt.subplot(2,2,2);
plt.imshow(T, interpolation = 'none', aspect = 'auto', cmap = 'viridis');
plt.colorbar(pad = 0.01, fraction = 0.01);
plt.title('temporal principal components')
plt.subplot(2,2,3);
plt.plot(100-pvar);
plt.title('variance explained')
plt.subplot(2,2,4, sharex = ax);
plt.xlim(0, len(dat_mean))
plt.plot(dat_mean)
plt.plot(dat_var)


reload(als)
reload(fplt)

plt.figure(10); plt.clf()
fplt.plot_pca(dist.T)

plt.figure(11); plt.clf();
fplt.plot_nmf(dist.T, n_components=None)

reload(fplt)
plt.figure(20);
fplt.plot_manifold_embeddings(dist.T, n_components=3, n_neighbors=20)
plt.tight_layout()


### Fit distributions

#dist_names = ['gamma', 'beta', 'rayleigh', 'norm', 'pareto'
import scipy;
dist_name = 'beta';

dist_form = getattr(scipy.stats, dist_name)

param = np.zeros((dist_form.numargs + 2, nbins));
eps = 10e-9;
for b in range(1):
  dd = dat_bin[:,b].copy();
  dd[dd > 0.7] = 0.7;
  dd = dd + eps;
  dd /= (dd.max() + eps);
  param[:, b] = dist_form.fit(dd, floc=0, fscale=1);

pdf_fitted = dist_form.pdf(np.linspace(1.0/dbins*0.5, 1-1.0/dbins, dbins), *param[:-2,b], loc=param[-2,b], scale=param[-1,b]);
pdf_fitted /= pdf_fitted.sum();

plt.figure(7); plt.clf();  
plt.plot(dist[:,b])
plt.plot(pdf_fitted, label=dist_name)


### Distribution for all scales

nscales = dat_scales.shape[1];
ntimes = dat_scales.shape[2];

dist_scales = np.zeros((nscales, dbins, ntimes));
dat_max_scales = np.zeros(nscales);
dat_min_scales = np.zeros(nscales);
for s in range(nscales):
  dat_max_scales[s] = np.percentile(dat_scales[:,s,:], 100);
  dat_min_scales[s] = np.percentile(dat_scales[:,s,:], 5);
  for t in range(ntimes):
    dist_scales[s,:,t] = np.histogram(dat_scales[:,s,t], range = (dat_min_scales[s],dat_max_scales[s]), bins = dbins)[0];
    dist_scales[s,:,t] /= dist_scales[s,:,t].sum();


dp = np.zeros((nscales*dbins, ntimes));
for i in range(nscales):
    dp[(i*dbins):((i+1)*dbins),:] = dist_scales[i,:,:]; 
cutoff = 0.2;
dp[dp > cutoff] = cutoff;
fplt.plot_array(dp)

plt.figure(10); plt.clf();
plt.plot(dat_max_scales,'b');
plt.plot(dat_min_scales,'r');





#plt.figure(11); plt.clf();
#fplt.plot_image_array(dist_scales[:5], names = [str(i) for i in range(5)])