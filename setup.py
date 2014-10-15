#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
# Using setuptools enables the package to be used in develop mode
# See http://pythonhosted.org/setuptools/setuptools.html#development-mode
# See https://github.com/scikit-learn/scikit-learn/issues/1016
try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    print('Please install or upgrade setuptools or pip to continue')
    sys.exit(1)


# See http://stackoverflow.com/a/18034855/2823213
have_cython = False
try:
    from Cython.Distutils import build_ext
    have_cython = True
except ImportError:
    from setuptools.command.build_ext import build_ext

if have_cython:
    _sim  = Extension('jittermodel._sim', ['jittermodel/_sim.pyx'])
else:
    _sim = Extension('jittermodel._sim', ['jittermodel/_sim.c'])

requirements = ['numpy', 'scipy', 'matplotlib', 'pint']

test_requirements = ['mpmath', 'bunch', 'nose']

setup(name='jittermodel',
      version='0.1',
      description='Calculate jitter and non-contact friction for an AFM cantilever',
      author='Ryan Dwyer',
      license='MIT',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7'],
      author_email='ryanpdwyer@gmail.com',
      url='http://github.com/ryanpdwyer/jittermodel',
      packages=['jittermodel', 'jittermodel.tests'],
      zip_safe=False,
      install_requires=requirements,
      tests_require=test_requirements,
      extras_require={
      'dev': ['sphinx']
      },
      ext_modules=[_sim],
      cmdclass={'build_ext': build_ext}
      )
