#!/usr/bin/env python3
# encoding: utf-8

from setuptools import setup

setup(
    name='XenoMapper',
    version='0.5.1b3',
    author='Matthew Wakefield',
    author_email='matthew.wakefield@unimelb.edu.au',
    packages=['xenomapper'],
    include_package_data = True,
    url='https://git@bitbucket.org/genomematt/xenomapper.git',
    license='GPL',
    entry_points={
        'console_scripts': ['xenomapper = xenomapper.xenomapper:main',
                            'xenomappability = xenomapper.mappability:main'
                           ]
    },
    test_suite='xenomapper.tests.test_all',
    description='xenomapper - mapping mixed reads from two species',
    long_description=open('README.md').read(),
    classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

)
