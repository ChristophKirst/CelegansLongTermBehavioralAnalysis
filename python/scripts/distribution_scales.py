# -*- coding: utf-8 -*-
"""
Sliding window distributions and scales

@author: ckirst
"""


import numpy as  np
import matplotlib.pyplot as plt;

import scipy.stats as stats

import experiment as exp
import analysis as als
import plot as fplt;


strain = 'cat-2';
feat = 'roam';

# calculate roaming a and dwelling intervalls -> ideally after smoothing with HHM ? 

data = exp.load_data(strain);

dat = getattr(data, feat);

#dat_scales = exp.load('%s_%s_time_normalized_scales_mean.npy' % (strain , feat));


dbins = 7;
sbins = 50;
nbins = 5 * sbins;

sbins = exp.stage_bins(data, nbins= sbins);

dat_bin = exp.bin_data(dat, sbins);
dat_mean = np.mean(dat_bin, axis = 0);
dat_var  = np.var(dat_bin, axis = 0);

dist = np.zeros((dbins, nbins));
for b in range(nbins):
  dist[:,b] += np.histogram(dat_bin[:,b], range = (0, 0.9), bins = dbins)[0];
  dist[:,b] /= np.sum(dist[:,b]);


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

def d(p,q):
  m = (p+q)/2;
  return stats.entropy(p,m) + stats.entropy(q,m);

ent2 = np.zeros((nbins, nbins));
for b in range(nbins):
  for b2 in range(nbins):
    ent2[b,b2] = d(dist[:,b], dist[:,b2]);
  


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