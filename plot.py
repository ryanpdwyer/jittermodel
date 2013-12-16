#!/usr/bin/env python
# encoding: utf-8
"""
plot.py

Created by Ryan Dwyer on 2013-10-15.
Copyright (c) 2013 Cornell University. All rights reserved.
"""
import numpy as np
from jittermodel.simulation import Simulation
from autoassign import autoassign
from time import sleep
import datetime
import cPickle as pickle
import os
import errno
import glob
import matplotlib.pyplot as plt


class GeneratePlotData(object):
    """An object that will generate plots of friction varying
    some variable var1 over the range var1_range.

    To initialize the object, pass a cantilever, sample and
    experiment object, variable and variable range to the object.
    Example:

    c1 = Cantilever()
    s1 = Sample()
    e1 = Experiment()
    fp1 = GenerateFrictionPlot(c1,s1,e1, 'V_g', )"""
    @autoassign
    def __init__(self, Cant, Samp, Expt, variable, variable_range):
        self.labels = {
            'rho': r'Carrier Density $\rho$ ($\mathrm{m}^{-3}$)',
            'V_g': r'Gate Voltage $V_\mathrm{g}$ (V)',
            'd': r'Tip-Sample Distance $d$ (nm)',
            'friction': (r'Friction $\Gamma_\mathrm{s}$ ' +
                         r'($\mathrm{pN} \mathrm{s}/\mathrm{m}$)'),
            'jitter': r'Jitter ($\mathrm{Hz}^2$)',
            'h': r'Sample Thickness $h$ (nm)',
            'power spectrum': (r'$P^{\perp}_{\delta f_c}(f)$ ' +
                               r'($\mathrm{Hz}^2 / \mathrm{Hz}$)'),
            'f': r'$f$ (Hz)'
        }
        self.scales = {
            'rho': 1e18,
            'V_g': 1e-3,
            'd': 1e3,
            'jitter': 1e6,
            'h': 1e3,
            'friction': 1,
            'power spectrum': 1e3,
            'f': 1e3
        }

    def _calc_variable_array(self, x_scale='log', n_pts=None):
        """Use the scale data to pick the appropriate
        spacing of the points."""
        if n_pts is None:
            n_pts = self.n_pts
        if x_scale == 'log':
            _r = np.logspace(np.log10(self.variable_range[0]),
                             np.log10(self.variable_range[1]), n_pts)
        elif x_scale == 'linear':
            _r = np.linspace(self.variable_range[0],
                             self.variable_range[1], n_pts)
        else:
            raise ValueError("x_scale must be either 'log' or 'linear'")
        return _r

    # def _pick_y_func(self, sim, label):
    #     """Picks the appropriate function from the class Simulation to calculate the labeled quantity. Currently, it accepts labels,
    #     
    #     friction
    #         Calculate the sample-induced friction using function calc_gamma_s
    #     
    #     jitter
    #         Calcualte the jitter
    #     """
    #     _func_dict = {'friction':sim.calc_gamma_s, 'jitter':sim.calc_jitter}
    #     try: 
    #         _func = _func_dict[label]
    #     except KeyError:
    #         raise KeyError('Please enter a correct label. The available labels are %s' % _func_dict.keys())
    #     return _func

    def _pick_plot_func(self, x_scale, y_scale):
        """Return the plot function with the appropriate scale on each axis."""
        _scale_tup = (x_scale, y_scale)
        _scale_dict = {('log', 'log'): plt.loglog,
                       ('linear', 'linear'): plt.plot,
                       ('linear', 'log'): plt.semilogy,
                       ('log', 'linear'): plt.semilogx}
        return _scale_dict[_scale_tup]

    # def plot(self, x_var, x_scale = 'log', y_scale = 'log', n_pts = 50):
    #     """This function plots the variable x_var, using the scales x_scale, y_scale, and using the number of points n_pts. The points at which x_var is evaluated are given by var1_range and var1, which are selected when the object is initialized."""
    #     
    #     """Convert the scale to a tuple, and then use a dictionary to pick the appropriate function for plotting."""
    #     _key_tup = (x_var, n_pts)
    #     
    #     var_range = self._calc_variable_array(x_scale, n_pts)
    #     
    #     """Takes the simulation data and varies var1 (defined when the object is initialized) over the range specified."""
    #     all_sims = []
    #     for val in var_range:
    #         sim = Simulation(self.Cant, self.Samp, self.Expt)
    #         sim[self.var1] = val
    #         # print sim
    #         all_sims.append(sim)
    #     
    #     
    #     """Saves the data to the object so that it is available for further use."""
    #     _func = self._pick_y_func(sim)
    #     self.plot_x[_key_tup] = _x_data = [sim[x_var] for sim in all_sims]
    #     self.plot_y[_key_tup] = _y_data = [sim.calc_gamma_s() for sim in all_sims]
    #     
    #     
    #     """Plot the results"""
    #     # Use LaTeX to render the axes labels.
    #     plt.rc('text', usetex=True)
    #     plt.rc('font', family='serif')
    #     _plot_func = self._pick_plot_func(x_scale, y_scale)
    #     _plot_func(_x_data, _y_data, 'b-')
    #     plt.xlabel(self.label_dict[x_var], fontsize = 8)
    #     plt.ylabel(r'Friction $\Gamma_\mathrm{s}$ ($\mathrm{pN} \mathrm{s}/\mathrm{m}$)', fontsize = 8)
    #     plt.savefig('fig_fric.pdf')

    def _make_sim_array(self):
        """Makes an array of simulations over the variable varied
        in the experiment, and also for multi_plot_var."""
        variable_values = self._calc_variable_array(self.x_scale)
        self.all_sims = []
        for multi_plot_val in self.multi_plot_values:
            row = []
            for val in variable_values:
                sim = Simulation(self.Cant, self.Samp, self.Expt)
                sim.assign(self.multi_plot_var, multi_plot_val)
                sim.assign(self.variable, val)
                row.append(sim)
            self.all_sims.append(row)

    def _calculate_plot_points(self):
        x = []
        y = []
        for sim_row in self.all_sims:
            x_row = [sim.lookup(self.x_var) for sim in sim_row]
            y_row = [sim.func_dict[self.y_var]() for sim in sim_row]
            x.append(x_row)
            y.append(y_row)
        self.x = np.array(x) * self.scales[self.x_var]
        self.y = np.array(y) * self.scales[self.y_var]

    @autoassign
    def calc_plot_data(self, y_var, x_var,
                       multi_plot_var, multi_plot_values,
                       x_scale='log', n_pts=50):
        """This plots y_var vs. x_var, over the range 'range'
        of 'variable' selected when the object was initialzied.
        y_var should be 'friction' or 'jitter'.

        Usage:

        gpd = GeneratePlotData(c1,s1,e1, 'V_g', (1e3, 50e3))
        gpd.calc_plot_data('friction', 'V_g', 'mobility', [1e-4, 1, 1e4])
        """
        self._make_sim_array()
        self._calculate_plot_points()

    def make_plot(self, figname=None, figsize=(3, 3), fontsize=10,
                  linewidth=2, xlim=None, ylim=None, properties=None):
        import matplotlib as mpl
        mpl.rcParams.update({
            'figure.figsize': figsize,
            'text.usetex': True,
            'font.family': 'serif',
            'font.serif': 'Times',
            'font.size': fontsize})
        import matplotlib.pyplot as plt
        lines = plt.loglog(self.x.T, self.y.T, linewidth=linewidth)
        plt.xlabel(self.labels[self.x_var])
        plt.ylabel(self.labels[self.y_var])

        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

        def reformat_properties(properties):
            """Reformats properties, given in the form,

            {'color':('b', 'g'), 'ls':('-','--')},

            into a format which can be looped over to apply
            them to individual lines."""
            reformatted_properties = [{} for i in properties.values()[0]]
            for key, tup in properties.items():
                for d, val in zip(reformatted_properties, tup):
                    d[key] = val
            return reformatted_properties

        if properties is not None:
            reformatted_properties = reformat_properties(properties)
            for d, line in zip(reformatted_properties, lines):
                plt.setp(line, **d)

        plt.tight_layout()
        if figname is None:
            plt.show()
            sleep(20)
            plt.close()
        else:
            plt.savefig(figname)
            plt.close()

    # def multi_plot_jitter(self, x_var, multi_plot_var, multi_plot_var, figname, x_scale = 'log', y_scale = 'log', n_pts = 50):
    #     """This calculates the jitter, integrated from f_i to f_f, and then"""
    #     scale_tup = (x_scale, y_scale)
    #     all_sims = self._make_sim_array(multi_plot_var, multi_plot_var_values, x_scale, n_pts)
    # 
    #     _multi_x_data = []
    #     _multi_y_data = []
    #     for sim_row in all_sims:
    #         x_row = [sim[x_var] for sim in sim_row]
    #         y_row = [sim.calc_jitter() for sim in sim_row]
    #         _multi_x_data.append(x_row)
    #         _multi_y_data.append(y_row)
    # 
    #     self.mplot_x  = _all_x = np.array(_multi_x_data)
    #     self.mplot_y = _all_y = np.array(_multi_y_data)

    def multi_power_spectrum(self, multi_plot_var, var_values,
                             f_i=1e-5, f_f=1e6, n_pts=100):
        """Plots the power spectrum, varying multi_plot_var over
        the values var_values. The user can also select the frequency
        range, which goes from an initial value of f_i to a final
        value of f_f. The user can also select the number of points at
        which to evaluate the power spectrum.


        TODO: Implement Marohn's protocol; save the data to a pickle file,
        and then read back the pickled data from a file in order to make the
        plot."""
        self.x_var = 'f'
        self.y_var = 'power spectrum'
        all_sims = []
        for multi_plot_val in var_values:
            sim = Simulation(self.Cant, self.Samp, self.Expt)
            sim.assign(multi_plot_var, multi_plot_val)
            all_sims.append(sim)

        _power_spectra_x = []
        _power_spectra_y = []
        for sim in all_sims:
            _x, _y = sim.calc_power_spectrum(f_i=f_i, f_f=f_f)
            _power_spectra_x.append(_x)
            _power_spectra_y.append(_y)

        self.x = np.array(_power_spectra_x)*self.scales['f']
        self.y = np.array(_power_spectra_y)*self.scales['power spectrum']

    def pickle(self, name):
        """Pickles the object to a file in a subdirectory pkl."""
        today = datetime.date.today().isoformat()
        filename = 'pkl/{today} {name}.pkl'.format(today=today, name=name)
        make_sure_path_exists('pkl')
        output = open(filename, 'wb')
        pickle.dump(self, output)
        output.close()


def unpickle(name=None):
    """Unpickles an object from the subdirectory pkl, with the given name,
    and created today."""
    if name is None:
        filename = get_most_recent_pkl_file('pkl')
    else:
        today = datetime.date.today().isoformat()
        filename = 'pkl/{today} {name}.pkl'.format(today=today, name=name)
    with open(filename, 'rb') as pkl_file:
        return pickle.load(pkl_file)


def get_most_recent_pkl_file(folder):
    """Gets the most recent pickle file from the specified folder.
    Code from http://stackoverflow.com/q/168409/2823213"""
    search_dir = "{cwd}/{folder}/*.pkl".format(cwd=os.getcwd(), folder=folder)
    # remove anything from the list that is not a file (directories, symlinks)
    # thanks to J.F. Sebastion for pointing out that the requirement was a list
    # of files (presumably not including directories)
    files = filter(os.path.isfile, glob.glob(search_dir))
    files.sort(key=lambda x: os.path.getmtime(x))
    return files[-1]


def make_sure_path_exists(path):
    """From Stackoverflow. See http://stackoverflow.com/q/273192/2823213"""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise  # This looks weird, but is correct. Allows any other error
                   # to be passed to the calling function.
