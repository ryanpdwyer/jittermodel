#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import jittermodel as jm
from jittermodel import u
from jittermodel.base import Cantilever, Transistor, Experiment
from jittermodel.plot import GeneratePlotData
import sys


if __name__ == '__main__':
    if sys.argv[1] == '1':
        cant = Cantilever(f_c=61.3*u.kHz, k_c=2.5*u.N/u.m,
                  Q=17980*u.dimensionless)
        trans = Transistor(mobility=5e-7*u('cm^2/V/s'), E_s2 = -0.5)
        expt = Experiment()

        gpd = GeneratePlotData(cant, trans, expt, 'd', (50*u.nm, 500*u.nm), model=1)

        V_g_vals = np.array([1, 50]) * u.V
        gpd.calc_plot_data('friction', 'd', 'V_g', V_g_vals, x_scale='log', n_pts=100)
    if sys.argv[1] == '2':
        cant = Cantilever(f_c=61.3*u.kHz, k_c=2.5*u.N/u.m,
                  Q=17980*u.dimensionless)
        trans = Transistor(mobility=5e-7*u('cm^2/V/s'), E_s2 = -0.5)
        expt = Experiment()

        gpd = GeneratePlotData(cant, trans, expt, 'd', (50*u.nm, 500*u.nm), model=2)

        V_g_vals = np.array([1, 50]) * u.V
        gpd.calc_plot_data('friction', 'd', 'V_g', V_g_vals, x_scale='log', n_pts=100)
