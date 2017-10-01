# -*- coding: utf-8 -*-
"""
Analyse Hidden Markov Models

http://pomegranate.readthedocs.io/en/latest/hmm.html
"""

import os 
import numpy as np
import matplotlib.pyplot as plt

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

def seq_from_path(path):
  return " ".join( state.name for idx, state in path[1:] )


def idx_from_path(path):
  return [idx for idx, state in path[1:]]

# make a sequence 

seq = [np.array(np.random.rand(100) > 0.2, dtype = int)];

model = hmm(nstates = 2);

nstates = 2;
states = [pmg.DiscreteDistribution({0 : 0.5, 1: 0.5}) for i in range(nstates)];
trans = np.ones((nstates, nstates)) / nstates;
trans = np.random.rand(nstates, nstates);
for i in range(nstates):
  trans[i] = trans[i] / trans[i].sum();
model = pmg.HiddenMarkovModel().from_matrix(trans, states, np.ones(nstates)/nstates, np.zeros(nstates));

model.plot()


print model.fit(seq)

plt.figure(1); plt.clf();
model.plot()

logp, path = model.viterbi(seq[0]);
print idx_from_path(path);




### worm data

import experiment as exp

strain = 'N2';
feat = 'roam';

data = exp.load_data(strain);

d = getattr(data, feat);

# train HMM on sequences in initial intervall

tstart = 450000;
tend = tstart + 6000;

sdat = np.array(d[:,tstart:tend], dtype = int);
model = hmm(nstates = 2);
model.fit(sdat, max_iterations = 500, stop_threshold =  10e-7, n_jobs = 12)



#logp, path = model.viterbi(d[wid, tstart:tend]);
#print seq_from_path(path)
wid = 10;
plt.figure(1); plt.clf();
plt.subplot(2,1,1);
model.plot();
plt.subplot(2,1,2);
plt.plot(model.predict(d[wid, tstart:tend]), 'r');
plt.plot(-d[wid, tstart:tend]+2.1, 'b')
#plt.plot(d[wid, tstart:tend] - idx_from_path(path))


## see how the hmm parameter change over time

def model_to_param(model):
  emis = [];
  for s in model.states[:-2]:
    emis.append(s.distribution.values()[-1]);
  trans = model.dense_transition_matrix();
  return np.hstack([trans[0,0], trans[1,1], emis, trans[2, :2]]);


sbins_start, sbins_end = exp.stage_bins(data, nbins = 2**5)

nbins_total = sbins_start.shape[1];
#nparam = nstates * (nstates-1) + nstates; #transition parameter + emission prob for 1
nparam = 6;
nworms = d.shape[0]
probs = np.zeros((nparam, nbins_total));
for i in range(nbins_total):
  print '%d/%d' % (i, nbins_total);
  model = hmm(nstates = 2);
  dd = [d[:,sbins_start[j][i]:sbins_end[j][i]] for j in range(nworms)];
  model.fit(dd, max_iterations = 500, stop_threshold =  10e-5, n_jobs = 12)
  probs[:,i] = model_to_param(model);
probs[3,:] = 1- probs[3,:]

droam = exp.bin_data(d, (sbins_start, sbins_end));

plt.figure(2); plt.clf();
plt.subplot(4,1,1);
plt.imshow(probs[:2,:], aspect = 'auto', interpolation = 'none', vmin = 0.8)
plt.subplot(4,1,2)
plt.imshow(probs[2:4,:], aspect = 'auto', interpolation = 'none', vmin = 0.85)
plt.subplot(4,1,3)
plt.imshow(probs[4:,:], aspect = 'auto', interpolation = 'none')
plt.subplot(4,1,4);
plt.plot(np.mean(droam, axis = 0))
plt.xlim(0, droam.shape[1])
plt.tight_layout()

np.save('%s_hmm.npy' % strain, probs)


probs = np.load('%s_hmm.npy' % strain)


### Single worm 

wid = 0;

model = hmm(nstates=2);

model.fit([d[wid]]);


model.fit(d, n_jobs = 12);


pred = 1- np.array(model.predict(d[wid]));


plt.figure(3); plt.clf();
ax = plt.subplot(2,1,1);
plt.plot(d[wid]);
plt.ylim(0, 1.5)
plt.title('data')
plt.subplot(2,1,2, sharex = ax);
plt.plot(pred);
plt.ylim(0, 1.5)                                                       
plt.title('pred')

import plot as fplt;
fplt.plot_array(np.vstack([d[wid], pred, d[wid]-pred]))


np.abs( d[wid]-pred).sum()

### isi distribution of this

import analysis as als


disi = als.isi_onoff(d)





