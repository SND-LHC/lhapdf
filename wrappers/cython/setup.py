#! /usr/bin/env python

from distutils.core import setup
from glob import glob
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension('lhapdf',
              ['lhapdf.pyx'],
              include_dirs=['../../include'],
              library_dirs=['../../src/.libs'],
              language='C++',
              libraries=['stdc++', 'LHAPDF']) ]

setup(name = 'lhapdflib',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules)
