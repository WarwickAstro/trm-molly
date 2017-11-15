trm.molly
=========

trm.molly is a Python-based package to access astronomical spectra written by
the F77 program 'molly'. trm.molly allows one to read and write molly format
data, and to get to the header and data in molly spectra from within Python.
The approach is very much the minimum "read the data in, make it available",
with anything more fancy left to the user.

trm.molly works under Python3.XX and needs the third-party packages 'numpy'
and 'scipy' to be installed.

Installation
------------

Once you have cloned trm.molly from github, you will have a sub-directory
called trm-molly. 'cd' to that, and then install in the usual way i.e.

python setup.py install

[might need super-user privileges, depending on your setup]

or

python setup.py install --user

to install locally [might need to adjust your PYTHONPATH to pick it up]


Let me know of problems / bugs,


Tom Marsh
