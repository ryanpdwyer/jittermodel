#!/usr/bin/env python
# encoding: utf-8
"""
scratchwork.py

Created by Ryan Dwyer on 2013-10-15.
Copyright (c) 2013 Cornell University. All rights reserved.
"""

from jittermodel.base import Cantilever, Sample, Experiment
from jittermodel.ubase import UnitCantilever, UnitTransistor, UnitExperiment
from jittermodel.simulation import Simulation
from jittermodel.plot import GeneratePlotData
from autoassign import autoassign
import pickle
import os
import errno

def make_sure_path_exists(path):
    """From Stackoverflow. See http://stackoverflow.com/questions/273192/create-directory-if-it-doesnt-exist-for-file-write for more information."""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


c1 = Cantilever()
s1 = Sample()
e1 = Experiment()
sim1 = Simulation(c1,s1,e1)
gpd1 = GeneratePlotData(c1, s1, e1, 'd', (40e-3,500e-3))
v_g_vals = (1e3, 20e3, 40e3)
gpd1.calc_plot_data('jitter', 'd', 'V_g', v_g_vals, n_pts = 5)
gpd1.make_plot('new_jitter_h_70.pdf', xlim = (40, 500), ylim = (1e-10, 1e-3))

cu1 = UnitCantilever()
tu1 = UnitTransistor()
eu1 = UnitExperiment()


# 
# 
# s2 = Sample(h = 25e-3)
# gpd2 = GeneratePlotData(c1, s2, e1, 'd', (40e-3,500e-3))
# gpd2.calc_plot_data('jitter', 'd', 'V_g', v_g_vals, n_pts = 5)
# gpd2.make_plot('jitter_h_25.pdf', xlim = (40, 500), ylim = (1e-10, 1e-3))
# 
# s3 = Sample(h = 10e-3)
# gpd3 = GeneratePlotData(c1, s2, e1, 'd', (40e-3,500e-3))
# gpd3.calc_plot_data('jitter', 'd', 'V_g', v_g_vals, n_pts = 5)
# gpd3.make_plot('jitter_h_10.pdf', xlim = (40, 500), ylim = (1e-10, 1e-3))
# 
# s4 = Sample(h = 3e-3)
# gpd4 = GeneratePlotData(c1, s2, e1, 'd', (40e-3,500e-3))
# gpd4.calc_plot_data('jitter', 'd', 'V_g', v_g_vals, n_pts = 5)
# gpd4.make_plot('jitter_h_3.pdf', xlim = (40, 500), ylim = (1e-10, 1e-3))


class Axis(object):
    @autoassign
    def __init__(self, variable, range, points):
        pass

    @property
    def min(self):
        return self.range[0]

    @property
    def max(self):
        return self.range[1]

    def __str__(self):
        return """{self.variable} on the interval \
                  {self.range} {self.points} points""".format(self=self)

h_ax = Axis('h', (10, 100), 50)

"""To pickle"""
make_sure_path_exists('pkl')
output = open('pkl/axis.pkl', 'wb')
pickle.dump(h_ax, output)
output.close()

# """To open pickled file"""
# pkl_file = open('pkl/axis.pkl', 'rb')
# h_ax_pkl = pickle.load(pkl_file)
# pkl_file.close()
# print h_ax_pkl

# d = {'a':1, 'b':2}
# output = open('data.pkl', 'wb')
# pickle.dump(d, output)
# output.close()

# pkl_file = open('data.pkl', 'rb')
# dd = pickle.load(pkl_file)

# print dd
