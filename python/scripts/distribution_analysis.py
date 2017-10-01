# -*- coding: utf-8 -*-
"""
Anaylsis of distributions

"""

import numpy as  np
import matplotlib.pyplot as plt;

import experiment as exp
import analysis as als
import plot as fplt;


strain = 'N2';
feat = 'roam';

# calculate roaming and dwelling intervalls -> ideally after smoothing with HHM ? 

data = exp.load_data(strain);

dat = getattr(data, feat);

# filter with hhm

import pomegranate as pmg

def hmm(nstates = 2, bias = 0.1):
  def make_bias(i,s):
    if i == 0:
      return [bias, 1-bias][s];
    else:
      return [1-bias, bias][s];
  states = [pmg.State(pmg.DiscreteDistribution({0 : make_bias(i,0), 1: make_bias(i,1)}), name = 'S%d'%i) for i in range(nstates)];
  
  #trans = np.ones((nstates, nstates)) / nstates;
  trans = np.random.rand(nstates, nstates);
  for i in range(nstates):
    trans[i] = trans[i] / trans[i].sum();
  
  model = pmg.HiddenMarkovModel();
  model.add_states(states);
  for i in range(nstates):
    for j in range(nstates):
      model.add_transition(states[i], states[j], trans[i,j]);
    model.add_transition(model.start, states[i], 1.0 / nstates );
  model.bake();
  return model;


model = hmm(nstates=2);
model.fit(dat, n_jobs = 12);


dat_hmm = np.zeros(dat.shape);
nworms = dat_hmm.shape[0];
for w in range(nworms):
  dat_hmm[w] = 1- np.array(model.predict(dat[w]));




# calculate duration epchs of romaing vs. dewlling ? -> why not summarize in terms of HHM parameter for each worm for some time window ?


(times_up, times_dw, dur_up, dur_dw) = als.isi_onoff(dat);

#(times_up, times_dw, dur_up, dur_dw) = als.isi_onoff(dat_hmm);


# make distributions

def bin_durations_for_worm(t_up, t_dw, d_up, d_dw, bin_st, bin_ed):
  nbins = len(bin_st);
  durs_up = []; 
  durs_dw = [];
  for b in range(nbins):
    idx_up = np.logical_and(bin_st[b] < t_up + d_up, t_up < bin_ed[b]);
    idx_dw = np.logical_and(bin_st[b] < t_dw + d_dw, t_dw < bin_ed[b]);
    durs_up.append(d_up[idx_up]);
    durs_dw.append(d_dw[idx_dw]);
  return (durs_up, durs_dw);
  
def bin_durations(t_up, t_dw, d_up, d_dw, bin_st, bin_ed):
  nworms = len(t_up);
  durs_up = [];
  durs_dw = [];
  for w in range(nworms):
    du,dd = bin_durations_for_worm(t_up[w], t_dw[w], d_up[w], d_dw[w], bin_st[w], bin_ed[w]);
    durs_up.append(du); durs_dw.append(dd);
  return (durs_up, durs_dw);


sbin_st, sbin_ed = exp.stage_bins(data, nbins=10);

durs_up, durs_dw = bin_durations(times_up, times_dw, dur_up, dur_dw, sbin_st, sbin_ed);


## distributions for individual worms

# durs have shape (nworms, nbins, xx)

nworms = len(durs_up);
nbins = len(sbin_st[0]);

# number of periods
ndurs_up = np.zeros((nworms, nbins));
ndurs_dw = np.zeros((nworms, nbins));
for w in range(nworms):
  for b in range(nbins):
    ndurs_up[w,b] = len(durs_up[w][b]);
    ndurs_dw[w,b] = len(durs_dw[w][b]);

# min/max values for each parameter

dur_up_max = np.zeros((nworms, nbins));
dur_up_min = np.zeros((nworms, nbins));
dur_dw_max = np.zeros((nworms, nbins));
dur_dw_min = np.zeros((nworms, nbins));

dur_up_mean = np.zeros((nworms, nbins));
dur_dw_mean = np.zeros((nworms, nbins));

for w in range(nworms):
  for b in range(nbins):
    du = durs_up[w][b];
    if len(du) > 0:
      dur_up_max[w,b] = du.max();
      dur_up_min[w,b] = du.min();
      dur_up_mean[w,b] = du.mean();
    else:
      dur_up_max[w,b] = 1;
      dur_up_min[w,b] = 0;
      dur_up_mean[w,b] = 0;
    
    du = durs_dw[w][b];
    if len(du) > 0:
      dur_dw_max[w,b] = du.max();
      dur_dw_min[w,b] = du.min();
      dur_dw_mean[w,b] = du.mean();
    else:
      dur_dw_max[w,b] = 1;
      dur_dw_min[w,b] = 0;
      dur_dw_mean[w,b] = 0;



dord = exp.load('%s_%s_order.npy' % (strain, feat));


mm = dur_up_mean.max()
dur_up_mean[dur_up_mean == mm] = 0;

plt.figure(2); plt.clf();
dur_names = ['dur_up_max', 'dur_up_min', 'dur_dw_max', 'dur_dw_min', 'dur_up_mean', 'dur_dw_mean', 'ndurs_up'];
fplt.plot_image_array((dur_up_max, dur_up_min, dur_dw_max, dur_dw_min, dur_up_mean, dur_dw_mean, ndurs_up), order = dord, names = dur_names)
plt.tight_layout()
# make histograms for each worm

dbins = 10;

dist_dur_up = np.zeros((nworms, nbins, dbins));
dist_dur_dw = np.zeros((nworms, nbins, dbins));

for w in range(nworms):
  for b in range(nbins):
    dist_dur_up[w,b] = np.histogram(durs_up[w][b], range = (0, dur_up_max[w,b]), bins = dbins)[0];
    dist_dur_dw[w,b] = np.histogram(durs_dw[w][b], range = (0, dur_dw_max[w,b]), bins = dbins)[0];


cutoff = 5;
pd_up = als.distributions_to_array(dist_dur_up[dord], worms_first=True);
pd_up[pd_up > cutoff] = cutoff;

pd_dw = als.distributions_to_array(dist_dur_dw[dord], worms_first=True);
pd_dw[pd_dw > cutoff] = cutoff;

plt.figure(4); plt.clf();
fplt.plot_image_array((pd_up, pd_dw), names = ['dist_dur_up', 'dist_dur_dw'], invert_y = True)



# cumulative

dist_dur_up_full = np.zeros((dbins, nbins));
dist_dur_dw_full = np.zeros((dbins, nbins));

dur_up_max_full = np.max(dur_up_max, axis = 0);
dur_dw_max_full = np.max(dur_dw_max, axis = 0);

dur_up_max_full = np.percentile(dur_up_max, 90, axis = 0);
dur_dw_max_full = np.percentile(dur_dw_max, 90, axis = 0);

dur_up_max_full = np.percentile(dur_up_max, 90, axis = 0);
dur_dw_max_full = np.percentile(dur_dw_max, 90, axis = 0);

for w in range(nworms):
  for b in range(nbins):
    #dist_dur_up_full[:,b] += np.histogram(durs_up[w][b], range = (0, dur_up_max_full[b]), bins = dbins)[0];
    #dist_dur_dw_full[:,b] += np.histogram(durs_dw[w][b], range = (0, dur_dw_max_full[b]), bins = dbins)[0];
    dist_dur_up_full[:,b] += np.histogram(durs_up[w][b], range = (0, dur_up_max_full[b]), bins = dbins)[0];
    dist_dur_dw_full[:,b] += np.histogram(durs_dw[w][b], range = (0, dur_dw_max_full[b]), bins = dbins)[0];

for b in range(nbins):
  dist_dur_up_full[:,b] /= np.sum(dist_dur_up_full[:,b]);
  dist_dur_dw_full[:,b] /= np.sum(dist_dur_dw_full[:,b]);

sbins = (sbin_st, sbin_ed);
dat_bin = exp.bin_data(dat, sbins);
dat_mean = np.mean(dat_bin, axis = 0);

plt.figure(5); plt.clf();
fplt.plot_image_array((dist_dur_up_full, dist_dur_dw_full, dat_bin), names = ['dist_dur_up_full', 'dist_dur_dw_full', 'mean'], invert_y = True, nplots = 6);
plt.subplot(6,1,4);
plt.plot(dur_up_max_full, 'r');
plt.title('dur_up_max_full')
plt.subplot(6,1,5)
plt.plot(dur_dw_max_full, 'b');
plt.title('dur_dw_max_full')
plt.subplot(6,1,6);
plt.plot(dat_mean)
plt.title('mean %s' % feat)





