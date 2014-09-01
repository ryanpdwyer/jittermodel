# encoding: utf-8
from jittermodel.plot import reformat_properties
from nose.tools import eq_


def test_reformat_properties():
    properties = {'color': ('b', 'g'), 'ls': ('-', '--')}
    exp_reformatted = [{'color': 'b', 'ls': '-'}, {'color': 'g', 'ls': '--'}]
    eq_(exp_reformatted, reformat_properties(properties))
