# -*- coding: utf-8 -*-
"""
scratchwork.py

Created by Ryan Dwyer on 2013-10-15.
Copyright (c) 2013 Cornell University. All rights reserved.
"""

from jittermodel import u
from jittermodel.base import Cantilever, Transistor, Experiment
from jittermodel.plot import GeneratePlotData
import cProfile
import pstats


def make_unit_plot():
    c1 = Cantilever()
    t1 = Transistor(h=70*u.nm)
    e1 = Experiment()
    gpd = GeneratePlotData(c1, t1, e1, 'd', (40*u.nm, 500*u.nm))
    V_g_vals = (1*u.V, 20*u.V, 40*u.V)
    gpd.calc_plot_data('jitter', 'd', 'V_g', V_g_vals,
                       x_scale='linear', n_pts=5)

cProfile.run('make_unit_plot()', 'unit_stats')
p = pstats.Stats('unit_stats')
