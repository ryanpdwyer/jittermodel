#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import subprocess
# Using setuptools enables the package to be used in develop mode
# See http://pythonhosted.org/setuptools/setuptools.html#development-mode
# See https://github.com/scikit-learn/scikit-learn/issues/1016
try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    print('Please install or upgrade setuptools or pip to continue')
    sys.exit(1)

# ------------------------------------------------------------------------------
# This section taken from
# http://nfoti.github.io/a-creative-blog-name/posts/2013/02/07/cleaning-cython-build-files/

args = sys.argv[1:]

# Make a `cleanall` rule to get rid of intermediate and library files
if "cleanall" in args:
    print "Deleting cython files..."
    # Just in case the build directory was created by accident,
    # note that shell=True should be OK here because the command is constant.
    subprocess.Popen("rm -rf build", shell=True, executable="/bin/bash")
    subprocess.Popen("rm -rf jittermodel/*.c", shell=True, executable="/bin/bash")
    subprocess.Popen("rm -rf jittermodel/*.so", shell=True, executable="/bin/bash")

    # Now do a normal clean
    sys.argv[1] = "clean"

# We want to always use build_ext --inplace
if args.count("build_ext") > 0 and args.count("--inplace") == 0:
    sys.argv.insert(sys.argv.index("build_ext")+1, "--inplace")

# ------------------------------------------------------------------------------

# See http://stackoverflow.com/a/18034855/2823213
have_cython = False
try:
    from Cython.Distutils import build_ext
    have_cython = True
except ImportError:
    from setuptools.command.build_ext import build_ext


def make_extension(name, have_cython):
    if have_cython:
        ext = '.pyx'
    else:
        ext = '.c'
    return Extension(name, [name.replace('.', '/')+ext])


extensions = ['jittermodel._sim']
ext_modules = [make_extension(name, have_cython) for name in extensions]

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
      ext_modules=ext_modules,
      cmdclass={'build_ext': build_ext},
      test_suite='nose.collector'
      )
