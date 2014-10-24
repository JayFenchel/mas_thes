#!/usr/bin/python
"""muAO-MPC: microcontroller applications online model predictive control.

Model Predictive Control (MPC) solver based on an augmented 
Lagrangian method in combination with a fast gradient
method. Input and state constraints are considered.
Automatic code generation with focus on embedded systems, 
with particular focus on microcontroller applications.
"""
# this file is based on numpy's setup.py, no copyright infringement intended. 
DOCLINES = __doc__.split("\n")

import os
from distutils.core import setup, Extension

import numpy.distutils.misc_util as np

import muaompc._ltidt.writecfiles as wf
import muaompc._ltidt as _ltidt

CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: Linux
Operating System :: MacOS
"""

NAME                = 'muAO-MPC'
MAINTAINER          = "Pablo Zometa"
MAINTAINER_EMAIL    = "pablo.zometa@ovgu.de"
DESCRIPTION         = DOCLINES[0]
LONG_DESCRIPTION    = "\n".join(DOCLINES[2:])
URL                 = "http://ifatwww.et.uni-magdeburg.de/syst/mpctool/"
LICENSE             = 'BSD'
CLASSIFIERS         = filter(None, CLASSIFIERS.split('\n'))
AUTHOR              = "Pablo Zometa, et.al."
AUTHOR_EMAIL        = "pablo.zometa@ovgu.de"
PLATFORMS           = ["Windows", "Linux", "Mac OS-X"]
MAJOR               = 0
MINOR               = 3
MICRO               = 0
ISRELEASED          = True
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

def write_version_py(filename='muaompc/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM MUAO-MPC SETUP.PY
version = '%(version)s'
"""
    with open(filename, 'w') as f:
        f.write(cnt % {'version': VERSION})
        
def get_package_data():
    # this is easier done using setuptools or the like, but for compatibility
    # we use only the standard library's distutils
    ltidt_base = 'muaompc/_ltidt'
    ltidt_data_dir = 'template'
    p_dir = {'muaompc._ltidt': ltidt_base}
    ltidt_data = []
    for dirname, dirnames, filenames in os.walk(ltidt_base + 
                                                '/' + ltidt_data_dir):
        for filename in filenames:
            rel_dirname = dirname.replace(ltidt_base + '/', '')
            ltidt_data.append(rel_dirname + '/' + filename)

    p_data = {'muaompc._ltidt': ltidt_data}
    return p_dir, p_data

def get_data_files():
    basedir = 'examples'
    f_data = []
    for dirname, dirnames, filenames in os.walk(basedir):
        for filename in filenames:
            f_data.append(dirname + '/' + filename)
    return basedir, f_data
    
def setup_package():
    wf._write_vla_c_files(_ltidt.__path__[0])
    write_version_py()
    
    _cmpc = Extension(
      'muaompc._cmpc',
      [
        'cmpc/mtx_ops.c',
        'cmpc/mpc_stc.c',
        'cmpc/mpc_ref.c',
        'cmpc/mpc_inc.c',
        'cmpc/mpc.c',
        'src/cmpc.c',
        'src/module.c'
      ],
      include_dirs = [
        'src',
        'cmpc/include'
      ] + np.get_numpy_include_dirs()
    )
    
    p_dir, p_data = get_package_data()
    f_dir, f_data = get_data_files()

    setup(name=NAME,
          version=VERSION,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          url=URL,
          license=LICENSE,
          classifiers=CLASSIFIERS,
          author=AUTHOR,
          author_email=AUTHOR_EMAIL,
          platforms=PLATFORMS,
          # test_suite='tests',
          ext_modules=[_cmpc],
          packages=['muaompc', 'muaompc._ltidt', 'muaompc._ltidt.tutorial'],
          package_dir = p_dir,
          package_data = p_data,
          data_files = [(f_dir, f_data), ('.', ['LICENSE.txt']), 
              ('muaompc/_ltidt/tutorial', ['muaompc/_ltidt/tutorial/Makefile',
                  'muaompc/_ltidt/tutorial/main_motor.c', 
                  'muaompc/_ltidt/tutorial/main_aircraft.c'])],
          )

if __name__ ==  '__main__':
    setup_package()
