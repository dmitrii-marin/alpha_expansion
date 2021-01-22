#File: setup.py
#!/usr/bin/python
from distutils.core import setup, Extension

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

pht_module = Extension('_alphaexpansion',
                       sources=['alpha_expansion_wrap.cxx',
                                'alpha_expansion.cpp',
                                'gco-v3.0/graph.cpp',
                                'gco-v3.0/maxflow.cpp',
                                'gco-v3.0/GCoptimization.cpp',
                                'gco-v3.0/LinkedBlockList.cpp',
                                ],
                       extra_compile_args = [
                          "-fopenmp",
                          "-std=c++11",
                          "-DGCO_ENERGYTYPE=double",
                       ],
                       include_dirs = [numpy_include, "."],
                       extra_link_args=['-lgomp']
                      )

setup(name = 'alphaexpansion',
        version = '0.1',
        author = 'Dmitrii Marin',
        description = 'Python wrapper for alpha-expansion algorithm',
        ext_modules = [pht_module],
        py_modules = ['alphaexpansion'],
    )
