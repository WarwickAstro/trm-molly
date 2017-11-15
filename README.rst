trm.molly
=========

trm.molly is a Python-based package to access astronomical spectra written by
the F77 program 'molly'. trm.molly allows one to read and write molly format
data, and to get to the header and data in molly spectra from within Python.
The approach is a minimal "read the data in, make it available", with anything
more fancy left to the user.

trm.molly needs the third-party packages 'numpy' and 'scipy' to be installed.
It works under Python3.XX and, I believe, Python2.7, but I have only tested it
under Python3 myself.

Installation
------------

Once you have cloned trm.molly from github, you will have a sub-directory
called trm-molly. 'cd' to that, and then install in the usual way i.e.

python setup.py install

[you might need super-user privileges, depending on your setup]

or

python setup.py install --user

to install locally [you might need to adjust your PYTHONPATH to pick it up]


Once you have installed it, 'pydoc trm.molly' should get you started.

Let me know of problems / bugs,


Tom Marsh
