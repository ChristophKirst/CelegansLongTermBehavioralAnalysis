# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 01:22:12 2016

@author: ckirst
"""

import experiment as exp
import plot as fplt

data = exp.load_data('N2');

data = exp.add_positions(data);


wid= 0;
fplt.plot_trace(data.positions[0], data.roam[0])