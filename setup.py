#!/usr/bin/env python
# Enables the package to be used in develop mode
# See http://pythonhosted.org/setuptools/setuptools.html#development-mode
# See https://github.com/scikit-learn/scikit-learn/issues/1016
try:
    import setuptools
except ImportError:
    pass

from distutils.core import setup

requirements = ['numpy', 'scipy', 'matplotlib']

test_requirements = ['mpmath', 'bunch', 'nosetests']

setup(name='jittermodel',
      version='0.1',
      description='',
      author='Ryan Dwyer',
      author_email='ryanpdwyer@gmail.com',
      packages=['jittermodel'],
      zip_safe=False,
      install_requires=requirements,
      tests_require=test_requirements
      )
