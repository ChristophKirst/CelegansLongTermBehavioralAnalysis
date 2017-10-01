# -*- coding: utf-8 -*-
"""
Experimental Data Access for Roaming and Dwelling Data

10/3/2016
"""
__license__ = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__author__ = 'Christoph Kirst <ckirst@rockefeller.edu>'
__docformat__ = 'rest'

import os
import scipy.io as io
import dill as pickle
import glob


import numpy as np
import matplotlib.pyplot as plt


strains = ['N2', 'cat-2', 'tdc-1', 'npr-1', 'daf-7', 'tph-1', 'ser-1', 'ser-4', 'ser-5', 'ser-7', 'mod-1', 'mod-5'];

dir_root = '/home/ckirst/Science/Projects/CElegans';

dir_base     = os.path.join(dir_root, 'Analysis/Roaming');
dir_analysis = os.path.join(dir_base, 'Data');
dir_fig      = os.path.join(dir_base, 'Figures');

dir_raw  = os.path.join(dir_root, 'Experiment/RawData_Roaming');

dir_data = os.path.join(dir_root, 'Experiment/Data_Roaming');
dir_pos  = os.path.join(dir_root, 'Experiment/Data_XY');


def load_raw(strain = 'N2', verbose = False):
  """Load experimental data from original mat files
  
    #column 1 - speed
    #column 2 - angular velocity 
    #column 5 - life stage
    #column 6 - roaming/dwelling classification (1- roaming 0-dwelling)
  """
  
  if strain == 'N2':
    # N2 has two data sets
    #mat = io.loadmat(os.path.join(dir_raw, 'N2_1.mat'));
    #data = mat['individual_Speed_AV2'][0];
    #mat2 = io.loadmat(os.path.join(dir_raw,'N2_2.mat'));
    #data2 = mat2['individual_Speed_AV2'][0];
    #data = np.hstack([data, data2]);
    #data2 = [];
    #mat = io.loadmat(os.path.join(dir_raw, 'N2_3.mat'));
    #data = mat['individual_Speed_AV2'][0];
    # N2 new data set
    mat = io.loadmat(os.path.join(dir_raw, 'N2_1_new.mat'));
    data = mat['N2_1'][0];
    mat2 = io.loadmat(os.path.join(dir_raw,'N2_2_new.mat'));
    data2 = mat2['N2_2'][0];
    data = np.hstack([data, data2]);
    data2 = [];  
  else:
    mat = io.loadmat(os.path.join(dir_raw, '%s.mat' % strain));
    data = mat['individual_Speed_AV2'][0];
  
  #create   
  d = lambda x: None;
  d.strain = strain;
  d.data = data;
  d.nworms = data.shape[0];
  d.total_time = [data[i].shape[0] for i in range(d.nworms)]
  d.max_time = max(d.total_time)  
  d.nstages = 5;
  
  d.speed = np.zeros((d.nworms, d.max_time));
  d.rot = np.zeros((d.nworms, d.max_time));
  d.roam  = np.zeros((d.nworms, d.max_time));
  d.stage = np.zeros((d.nworms, d.max_time));
  for i in range(d.nworms):
    d.speed[i,0:d.total_time[i]] = data[i][:,0];
    d.rot[  i,0:d.total_time[i]] = data[i][:,1];
    d.roam[ i,0:d.total_time[i]] = data[i][:,5];
    d.stage[i,0:d.total_time[i]] = data[i][:,4];
  
  ### Stage Analysis
  d.stage_switch = np.zeros((d.nworms,d.nstages+1), dtype = int);
  for i in range(d.nworms):
    d.stage_switch[i,1:-1] = np.where(np.diff(data[i][:,4])>0)[0];
    d.stage_switch[i,-1] = d.total_time[i];
  
  d.stage_durations = np.diff(d.stage_switch, axis = 1);
  d.stage_min_durations = np.min(d.stage_durations, axis = 0);
  
  if verbose:
    plt.figure(1); plt.clf();
    for i in range(d.nstages):
      plt.hist(d.stage_durations[:,i], bins = 200, range = [0, d.stage_durations.max()]);
    plt.title('stage durations')
    
  return d;


def dataname(name):
  return os.path.join(dir_data, name);

def figname(name):
  return os.path.join(dir_fig, name);

def load_data(strain):
  """Load data from pickeled file"""
  return pickle.load(open(os.path.join(dir_data, '%s.pickle' % strain), 'rb'));
  
def save_data(data):
  """Save data to pickeled file"""
  pickle.dump(data, open(os.path.join(dir_data, '%s.pickle' % data.strain), 'wb'));


def load_analysis(name):
  """Load Data from analysis directory"""
  if len(name) > 4 and name[-4:] == '.npy':
    return np.load(os.path.join(dir_data, name), mmap_mode = memmap);
  else:
    return pickle.load(open(os.path.join(dir_data, name), 'rb'));

def save_analysis(data, name):
  """Save data to analysis directory"""
  if isinstance(data,np.ndarray):
    np.save(os.path.join(dir_data, name), data)  ;
  else:
    pickle.dump(data, open(os.path.join(dir_data, name), 'wb'));

  
def load(name, memmap = None):
  """Load Data from experiment directory"""
  if len(name) > 4 and name[-4:] == '.npy':
    return np.load(os.path.join(dir_data, name), mmap_mode = memmap);
  else:
    return pickle.load(open(os.path.join(dir_data, name), 'rb'));

def save(data, name):
  """Save data to experiment directory"""
  if isinstance(data,np.ndarray):
    np.save(os.path.join(dir_data, name), data)  ;
  else:
    pickle.dump(data, open(os.path.join(dir_data, name), 'wb'));


def ls(name, print_files = True, remove_directory = True):
  """Return list of the files that match name"""
  fns = glob.glob(os.path.join(dir_data, name));
  if remove_directory:
    fns = [os.path.split(f)[1] for f in fns];
  if print_files:
    print "\n".join(fns);
  else:
    return fns;


def memmap(name, shape, mode = 'r', dtype = None):
  """Create new memmap in experiment directory
  Note:
    This is usefull for parallel processing
  """
  filename = os.path.join(dir_data, name)
  return np.lib.format.open_memmap(filename, mode = mode, shape = shape);


def convert_data(strain):
  """Convert data from mat to pickle for faster access"""
  print 'converting strain %s' % strain;
  data = load_raw(strain);
  save_data(data);


def load_positions(strain, wid):
  fn = os.path.join(dir_pos, '%s_xy_w=%d_s=all.npy' % (strain.lower(), wid));
  return np.load(fn);

#def load_image(strain, wid, t):
#  fn = os.path.join(dir_pos, '%s_img_w=%d_s=all.npy' % (strain.lower(), wid));
#  fimg = np.load(fn, mmap_mode = 'r');
#  return fimg[t];

def add_positions(data):
  pos = np.zeros((data.nworms, data.max_time, 2));
  for wid in range(data.nworms):
    p = load_positions(data.strain, wid);
    pos[wid, :p.shape[0],:] = p;
  data.positions = pos;
  return data;


def stage_bins(data, nbins = 10, normalized = True):
  """Calcualte time normalized bins for the data set, nbins per stage"""
  nbins_total = nbins * data.nstages;
  
  if normalized:
    #normalized bins
    dt = 1.0 * data.stage_durations / nbins
    bins = np.zeros((data.nworms, nbins_total+1), dtype = int);
    for i in range(data.nworms):
      for s in range(data.nstages):
        bins[i, s*nbins:(s+1)*nbins] = np.asarray(np.arange(0,nbins) * dt[i,s] + data.stage_switch[i,s], dtype = int);
      bins[i,-1] = data.total_time[i];
    
    bins_start = bins[:,:-1];
    bins_end   = bins[:,1:];
    
    return bins_start, bins_end;
  
  else:
    #unormalized bins
    dt = 1.0 * data.stage_min_durations / nbins;
    ubins_start = np.zeros((data.nworms, nbins_total), dtype = int);
    ubins_end   = np.zeros((data.nworms, nbins_total), dtype = int);
    for i in range(data.nworms):
      for s in range(data.nstages):
        ubins_start[i, s*nbins:(s+1)*nbins] = np.asarray(np.arange(0,nbins) * dt[s] + data.stage_switch[i,s], dtype = int);
        ubins_end[i, s*nbins:(s+1)*nbins] = np.asarray(np.arange(1,nbins+1) * dt[s] + data.stage_switch[i,s], dtype = int);
    
    return ubins_start, ubins_end;
  

def bin_data(d, bins):
  """Bin the data according to the specified ranges"""
  bins_start, bins_end = bins;
  nbins_total = bins_start.shape[1];
  binned = np.zeros((d.shape[0], nbins_total));
  for i in range(d.shape[0]):
    binned[i,:] = np.array([np.mean(d[i,s:e]) for s,e in zip(bins_start[i], bins_end[i])]);
  return binned;


def align_stages(d, zeros):
  """Align data to specific stages"""
  nworms,ntimes = d.shape;
  lleft = max(zeros); 
  lright = max(ntimes - zeros);
  l = lleft + lright;
  a = np.zeros((nworms,l), dtype = d.dtype);
  for i in range(nworms):
    a[i,lleft-zeros[i]:lleft-zeros[i]+ntimes] = d[i,:];
  return a;



if __name__ == '__main__':
  
  import experiment as exp;
  for s in exp.strains:
    exp.convert_data(s);

