from setuptools import find_packages, setup

import os

from Cython.Build import cythonize

with open('README.md') as f:
    long_description = f.read()

os.environ['CFLAGS'] = '-std=c++2a'

setup(name="drawtetrado",
      version="1.2",
      packages=['drawtetrado'],
      package_dir={'': 'src'},
      author="Michal Zurkowski",
      author_email="michal.zurkowski@cs.put.poznan.pl",
      description="Draw simplified, layer diagrams of quadruplexes.",
      long_description=long_description,
      long_description_content_type='text/markdown',
      url="https://github.com/michal-zurkowski/drawtetrado",
      project_urls={
          'Bug Tracker': 'https://github.com/michal-zurkowski/drawtetrado/issues'
      },
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      ext_modules=cythonize(["cython/optimizer.pyx"]),
      entry_points={'console_scripts': ['drawtetrado=drawtetrado.main:main']},
      install_requires=['pycairo', 'svgwrite'])
