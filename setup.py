import os

from setuptools import setup
from Cython.Build import cythonize

os.environ['CFLAGS'] = '-std=c++20'
setup(ext_modules=cythonize("optimizer.pyx"))
