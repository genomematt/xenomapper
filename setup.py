#!/usr/bin/env python3
# encoding: utf-8

from setuptools import setup

install_requires = []

try:
    import statistics
except ImportError:
    install_requires.append('statistics')

setup(
    name='XenoMapper',
    version='1.0.0',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    packages=['xenomapper'],
    include_package_data = True,
    install_requires=install_requires,
    url='https://github.com/genomematt/xenomapper.git',
    license='GPLv3',
    entry_points={
        'console_scripts': ['xenomapper = xenomapper.xenomapper:main',
                            'xenomappability = xenomapper.mappability:main'
                           ]
    },
    test_suite='xenomapper.tests.test_all',
    description='xenomapper - mapping mixed reads from two species',
    long_description=open('README.md').read(),
    classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3 :: Only',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: Implementation :: CPython',
          'Programming Language :: Python :: Implementation :: PyPy',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

)
