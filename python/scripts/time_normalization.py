# -*- coding: utf-8 -*-
"""
Generate time normalized data sets
"""
import numpy as np
import matplotlib.pyplot as plt

import experiment as exp
import plot as fplt

strain = 'N2';
feat = 'roam';

nbins = 2**13;

save_fig = True;

data = exp.load_data(strain);



### Bin data

sbins = exp.stage_bins(data, nbins = nbins);

d = getattr(data, feat);
tn = exp.bin_data(d, sbins);

exp.save(tn, '%s_%s_time_normalized.npy' % (strain, feat))

#fig = plt.figure(1); plt.clf();
#plt.imshow(tn, aspect = 'auto');
#plt.tight_layout();
#plt.title('%s %s time normalized' % (strain, feat));
#if save_fig:
#  fig.savefig(exp.figname('%s_%s_time_normalized.png'% (strain, feat)));

fig = fplt.plot_array(tn, title = '%s %s time normalized' % (strain, feat));
fplt.savefig(fig,exp.figname('%s_%s_time_normalized.png'% (strain, feat)), width = 2500)