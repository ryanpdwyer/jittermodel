# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
# See http://stackoverflow.com/a/3054314
matplotlib.use('Agg')

from nose.tools import eq_
import unittest


from jittermodel import u, silentremove
from jittermodel.base import Cantilever, Transistor, Experiment
from jittermodel.plot import reformat_properties, GeneratePlotData, unpickle
from jittermodel.tests import expected_failure


class Test_reformat_properties(unittest.TestCase):

    @staticmethod
    def test_reformat_properties():
        properties = {'color': ('b', 'g'), 'ls': ('-', '--')}
        exp_reformatted = [{'color': 'b', 'ls': '-'},
                           {'color': 'g', 'ls': '--'}]
        eq_(exp_reformatted, reformat_properties(properties))


class TestGeneratePlotData(unittest.TestCase):
    filename = 'test-generate-plot-data.pkl'
    figname = 'test-generate-plot-data-plot.png'

    def setUp(self):
        # If I actually care what is executed, replace these with objects where
        # every possible parameter is specified.
        cant = Cantilever()
        trans = Transistor()
        expt = Experiment()
        self.gpd = GeneratePlotData(cant, trans, expt, 'd', (40*u.nm, 500*u.nm))

    def tearDown(self):
        silentremove(self.filename)
        silentremove(self.figname)

    def test_plot_linear_scale(self):
        # Setup the plotting object without actually calling the plotting
        # function, because we don't want to waste time actually calculating
        # friction values.
        self.gpd.x_var = 'd'
        self.gpd.y_var = 'friction'
        self.gpd.multi_plot_var = 'V_g'
        self.gpd.multi_plot_values = [1*u.V]
        self.gpd.x_scale = 'linear'
        self.gpd.n_pts = 50
        self.gpd.x = np.linspace(40, 500, 50)  # Generic x data
        self.gpd.y = self.gpd.x * 2  # Generate y data
        self.gpd.make_plot(figname=self.figname, tight_layout=False)

    def test_pickle_GeneratePlotData_with_no_plot_data(self):
        """The simulation has not been run yet, so the pickling should work
        fine."""
        self.gpd.pickle(self.filename)
        unpickle(self.filename)
