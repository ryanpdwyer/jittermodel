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

import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = 'jittermodel/_version.py'
versioneer.versionfile_build = 'jittermodel/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = 'jittermodel-' # dirname like 'myproject-1.2.0'


# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------

# This section taken from
# http://nfoti.github.io/a-creative-blog-name/posts/2013/02/07/cleaning-cython-build-files/

args = sys.argv[1:]

# Make a `cleanall` rule to get rid of intermediate and library files
if "cleanall" in args:
    print("Deleting cython files...")
    import glob
    import shutil
    import os.path
    try:
        shutil.rmtree('build')
    except:
        pass

    for filename in glob.glob(os.path.join('jittermodel','*.c')):
        os.remove(filename)

    for filename in glob.glob(os.path.join('jittermodel','*.so')):
        os.remove(filename)

    # Now do a normal clean
    sys.argv[1] = "clean"

# We want to always use build_ext --inplace
if "build_ext" in sys.argv and "--inplace" not in sys.argv:
    sys.argv.insert(sys.argv.index("build_ext")+1, "--inplace")

# See http://stackoverflow.com/a/18034855/2823213
# If the user does not have Cython installed, we should be able to use the
# package by building with the .c file
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

# List of cython modules here. Make_extension will add the correct extension.
# Assumes the file lives at pkgname/modulename
extensions = ['jittermodel._sim']
ext_modules = [make_extension(name, have_cython) for name in extensions]

# Versioneer commands
cmdclass=versioneer.get_cmdclass()
# Add Cython / setuptools build_ext
cmdclass.update({'build_ext': build_ext})

setup(name='jittermodel',
      version=versioneer.get_version(),
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
      packages=find_packages(),
      zip_safe=False,
      install_requires=['numpy', 'scipy', 'matplotlib', 'pint'],
      tests_require=['mpmath', 'bunch', 'nose'],
      extras_require={
      'dev': ['sphinx']
      },
      ext_modules=ext_modules,
      cmdclass=cmdclass,
      test_suite='nose.collector',
      )
