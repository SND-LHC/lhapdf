from distutils.core import setup
from glob import glob
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension('lhapdf',
              ['lhapdf.pyx'] + glob('../../src/*.cpp'),
              include_dirs=['../../include'],
              language='C++',
              libraries=['stdc++']) ]

setup(name = 'lhapdflib',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules)
