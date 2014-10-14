# -*- coding: utf-8 -*-
"""
plot.py

Created by Ryan Dwyer on 2013-10-15.
Copyright (c) 2013 Cornell University. All rights reserved.
"""
import numpy as np
from jittermodel.simulation import Simulation
from time import sleep
import datetime
import cPickle as pickle
import os
import errno
import glob
import matplotlib.pyplot as plt
from jittermodel import u


def reformat_properties(properties):
    """Reformats properties, given in the form,

    properties = {'color':('b', 'g'), 'ls':('-','--')}

    into a format which can be looped over to apply
    them to individual lines."""
    reformatted_properties = [{} for i in properties.values()[0]]
    for key, tup in properties.items():
        for d, val in zip(reformatted_properties, tup):
            d[key] = val
    return reformatted_properties


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
    Simulation = Simulation
    units = Simulation.units

    def __init__(self, Cant, Samp, Expt, variable, variable_range, model=2):
        self.Cant = Cant
        self.Samp = Samp
        self.Expt = Expt
        self.variable = variable
        self.variable_range = variable_range
        self.model = model
        self.labels = {
            'rho': r'Carrier Density $\rho$ ($\mathrm{cm}^{-3}$)',
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
        # Shouldn't scales be inferred from self.units, self.output units?
        # After looking at it, yes it should.
        self.scales = {
            'rho': 1e12,
            'V_g': 1e-3,
            'd': 1e3,
            'jitter': 1e6,
            'h': 1e3,
            'friction': 1,
            'power spectrum': 1e3,
            'f': 1e3
        }

        self.output_units = {
            'rho': u.cm**-3,
            'V_g': u.V,
            'd': u.nm,
            'jitter': u.Hz**2,
            'h': u.nm,
            'friction': u.pN*u.s/u.m,
            'power spectrum': u.Hz,
            'f': u.Hz}

    def _calc_variable_array(self, x_scale='log', n_pts=None):
        """Use the scale data to pick the appropriate
        spacing of the points."""
        x_unit = self.output_units[self.x_var]
        x_min = self.variable_range[0].to(x_unit).magnitude
        x_max = self.variable_range[1].to(x_unit).magnitude

        if n_pts is None:
            n_pts = self.n_pts
        if x_scale == 'log':
            x_unit = self.output_units[self.x_var]
            x_min_log = np.log10(x_min)
            x_max_log = np.log10(x_max)
            x = np.logspace(x_min_log, x_max_log, n_pts) * x_unit
        elif x_scale == 'linear':
            x = np.linspace(x_min, x_max, n_pts) * x_unit
        else:
            raise ValueError("x_scale must be either 'log' or 'linear'")
        return x

    def _pick_plot_func(self, x_scale, y_scale):
        """Return the plot function with the appropriate scale on each axis."""
        scales_tuple = (x_scale, y_scale)
        scales_plotting_function_dict = {('log', 'log'): plt.loglog,
                       ('linear', 'linear'): plt.plot,
                       ('linear', 'log'): plt.semilogy,
                       ('log', 'linear'): plt.semilogx}
        return scales_plotting_function_dict[scales_tuple]

    def _make_sim_array(self):
        """Makes an array of simulations over the variable varied
        in the experiment, and also for multi_plot_var."""
        x_array = self._calc_variable_array(self.x_scale)
        self.all_sims = []
        for multi_plot_val in self.multi_plot_values:
            row = []
            for x in x_array:
                sim = self.Simulation(self.Cant, self.Samp, self.Expt,
                                      self.model)
                sim.assign(self.multi_plot_var, multi_plot_val)
                sim.assign(self.variable, x)
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
        self.y_var = y_var
        self.x_var = x_var
        self.multi_plot_var = multi_plot_var
        self.multi_plot_values = multi_plot_values
        self.x_scale = x_scale
        self.n_pts = n_pts

        self._make_sim_array()
        self._calculate_plot_points()

    def make_plot(self, figname=None, figsize=(3, 3), fontsize=10,
                  linewidth=2, xlim=None, ylim=None, properties=None,
                  tight_layout=True, transparent=False, usetex=False):
        import matplotlib as mpl

        if usetex:
            mpl.rcParams.update({
                                'figure.figsize': figsize,
                                'text.usetex': True,
                                'font.family': 'serif',
                                'font.serif': 'Times',
                                'font.size': fontsize})
        else:
            mpl.rcParams.update({
                                'figure.figsize': figsize,
                                'font.size': fontsize})

        import matplotlib.pyplot as plt
        plot = self._pick_plot_func(self.x_scale, 'log')
        lines = plot(self.x.T, self.y.T, linewidth=linewidth)
        plt.xlabel(self.labels[self.x_var])
        plt.ylabel(self.labels[self.y_var])

        if xlim is not None:
            plt.xlim(xlim)
        if ylim is not None:
            plt.ylim(ylim)

        if properties is not None:
            reformatted_properties = reformat_properties(properties)
            for d, line in zip(reformatted_properties, lines):
                plt.setp(line, **d)

        if tight_layout:
            plt.tight_layout()

        if figname is None:
            plt.show()
            sleep(20)
            plt.close()
        else:
            plt.savefig(figname, transparent=transparent)
            plt.close()

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
            sim = self.Simulation(self.Cant, self.Samp, self.Expt)
            sim.assign(multi_plot_var, multi_plot_val)
            all_sims.append(sim)

        power_spectra_x = []
        power_spectra_y = []
        for sim in all_sims:
            x, y = sim.calc_power_spectrum(f_i=f_i, f_f=f_f)
            power_spectra_x.append(x)
            power_spectra_y.append(y)

        self.x = np.array(power_spectra_x)*self.scales['f']
        self.y = np.array(power_spectra_y)*self.scales['power spectrum']

    def pickle(self, name):
        """Pickles the object to a file in a subdirectory pkl."""
        self.Cant = self.Expt = self.Samp = self.Simulation = None
        today = datetime.date.today().isoformat()
        filename = 'pkl/{today}-{name}.pkl'.format(today=today, name=name)
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
        filename = 'pkl/{today}-{name}.pkl'.format(today=today, name=name)
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
