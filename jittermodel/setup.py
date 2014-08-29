try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Plot CV',
    'author': 'Ryan Dwyer',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'ryanpdwyer@gmail.com',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['plot_cv'],
    'scripts': [],
    'name': 'projectname'
}

setup(**config)