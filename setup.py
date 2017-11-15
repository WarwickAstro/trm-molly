from setuptools import setup

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(

    name='trm.molly',
    version='1',
    packages = ['trm','trm.molly',],

    description="Access to molly spectra from Python",
    long_description=long_description,

    # metadata
    author='Tom Marsh',
    author_email='t.r.marsh@warwick.ac.uk',
    url='http://www.astro.warwick.ac.uk/',

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy','scipy'],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts' : [
            'genspec=trm.molly.command_line:genspec',
        ],
    },

)

