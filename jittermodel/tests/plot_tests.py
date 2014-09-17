# encoding: utf-8
import numpy as np

from jittermodel import u, silentremove
from jittermodel.base import Cantilever, Transistor, Experiment
from jittermodel.plot import reformat_properties, GeneratePlotData, unpickle
from nose.tools import eq_


def test_reformat_properties():
    properties = {'color': ('b', 'g'), 'ls': ('-', '--')}
    exp_reformatted = [{'color': 'b', 'ls': '-'}, {'color': 'g', 'ls': '--'}]
    eq_(exp_reformatted, reformat_properties(properties))


# class TestGeneratePlotData():
#     filename = 'test-generate-plot-data.pkl'

#     def setUp(self):
#         cant = Cantilever()
#         trans = Transistor()
#         expt = Experiment()
#         self.gpd = GeneratePlotData(cant, trans, expt, 'd', (40*u.nm, 500*u.nm))
#         V_g_vals = np.array([1, 5, 25]) * u.V
#         self.gpd.calc_plot_data('friction', 'd', 'V_g', V_g_vals)

#     def tearDown(self):
#         silentremove(self.filename)

#     def test_pickle_GeneratePlotDataObject(self):
#         self.gpd.pickle(self.filename)
#         unpickle(self.filename)
