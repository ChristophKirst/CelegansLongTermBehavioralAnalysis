# -*- coding: utf-8 -*-
"""
Analysis tools for the Roaming/Dwelling data set
"""
__license__ = 'MIT License <http://www.opensource.org/licenses/mit-license.php>'
__author__ = 'Christoph Kirst <ckirst@rockefeller.edu>'
__docformat__ = 'rest'


import numpy as np

import scipy.stats as stats


def jensen_shannon_divergence(p,q):
  """Jensen-Shannon distance between distributions p and q"""
  m = (p+q)/2.0;
  return stats.entropy(p,m) + stats.entropy(q,m);

def bhattacharyya_distance(p,q):
  return np.sqrt(1 - np.sqrt(np.sum(p * q)))

def pca(Y):
  """PCA with temporal (T) and spatial eigenvectors (C)"""
  u,s,C = np.linalg.svd(np.dot(Y,Y.T));
  Winv = np.diag(1.0/np.sqrt(s));
  L = np.dot(Winv, u.T);
  T = np.dot(L, Y);
  pvar = 100*(s)/s.sum();
  return C,T,pvar
 
 
def scales_to_array(data, worms_first = False, order = None):
  """Convert time scale resolved data into image"""
  nworms, nscales, ntimes = data.shape;
  res = np.zeros((nworms*nscales, ntimes));
  if order is None:
    order = range(nworms);
  if worms_first:
    for i in range(nscales):
      res[(i*nworms):((i+1)*nworms),:] = data[order,i,:];
  else:
    for i in range(nworms):
      res[(i*nscales):((i+1)*nscales),:] = data[order[i],:,:];
  return res;
  
  
def distributions_to_array(data, worms_first = False):
  """Convert distributions for worms andtimes into image"""
  nworms, ntimes, ndist = data.shape;
  res = np.zeros((nworms*ndist, ntimes));
  if worms_first:
    for i in range(ndist):
      res[(i*nworms):((i+1)*nworms),:] = data[:,:,i]; 
  else:
    for i in range(nworms):
      res[(i*ndist):((i+1)*ndist),:] = data[i,:,:].T;
  return res;


def isi_onoff(data):
  """Calculate ISIs for a 0/1 classifications of worm time traces"""
  if data.ndim > 1:
    nworms,ntimes = data.shape;
  else:
    data = np.array([data]);
    nworms = 1;
  
  sw = np.diff(data, axis = 1);
  sw[:,0] = -1; # make a down transition initially
  
  dur_up = [];
  times_up = [];
  dur_dw = [];
  times_dw = [];
  for i in range(nworms):
    t_up = np.where(sw[i,:] == 1)[0];  # ups
    t_dw = np.where(sw[i,:] == -1)[0]; # downs
    if len(t_dw) > len(t_up):          #last intervall assumed to be cutoff and thus neglected
      dur_up.append(t_dw[1:] - t_up);
      dur_dw.append(t_up - t_dw[:-1]);
      times_up.append(t_up);
      times_dw.append(t_dw[:-1]);
    else:
      dur_up.append(t_dw[1:] - t_up[:-1]);
      dur_dw.append(t_up - t_dw);
      times_up.append(t_up[:-1]);
      times_dw.append(t_dw);
  
  return (times_up, times_dw, dur_up, dur_dw);
  
  
    
def correctPValues(pvalues, method = 'BH'):
    """Corrects p-values for multiple testing using various methods 
    
    Arguments:
        pvalues (array): list of p values to be corrected
        method (Optional[str]): method to use: BH = FDR = Benjamini-Hochberg, B = FWER = Bonferoni
    
    Returns:
        array: correctedf p-values (q-values)
    
    References: 
        - `Benjamini Hochberg, 1995 <http://www.jstor.org/stable/2346101?seq=1#page_scan_tab_contents>`_
        - `Bonferoni correction <http://www.tandfonline.com/doi/abs/10.1080/01621459.1961.10482090#.VmHWUHbH6KE>`_
        - `R statistics package <https://www.r-project.org/>`_
    
    Notes:
        - modified from http://statsmodels.sourceforge.net/ipdirective/generated/scikits.statsmodels.sandbox.stats.multicomp.multipletests.html
    """
    
    pvals = np.asarray(pvalues);

    if method.lower() in ['bh', 'fdr']:
        
        pvals_sorted_ids = np.argsort(pvals);
        pvals_sorted = pvals[pvals_sorted_ids]
        sorted_ids_inv = pvals_sorted_ids.argsort()

        n = len(pvals);
        bhfactor = np.arange(1,n+1)/float(n);

        pvals_corrected_raw = pvals_sorted / bhfactor;
        pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
        pvals_corrected[pvals_corrected>1] = 1;
    
        return pvals_corrected[sorted_ids_inv];
    
    elif method.lower() in ['b', 'fwer']:
        
        n = len(pvals);        
        
        pvals_corrected = n * pvals;
        pvals_corrected[pvals_corrected>1] = 1;\
        
        return pvals_corrected;
        
    #return reject[pvals_sortind.argsort()]