__author__ = 'jayf'
from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Own math functions',
  ext_modules = cythonize("MyMath/MyMath.pyx"),
)