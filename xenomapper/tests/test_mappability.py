#!/usr/bin/env python3
# encoding: utf-8
"""
test_mappability.py

Created by Matthew Wakefield on 2012-08-16.
Copyright (c) 2012 Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""

import unittest
import sys, io
from xenomapper.mappability import simulate_reads
import hashlib
from pkg_resources import resource_stream

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2015 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.5.1b3"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development"

class test_mappability(unittest.TestCase):
    def setUp(self):
        pass
    def test_simulate_reads(self):
        outputfasta = io.StringIO()
        inputfastafile = io.TextIOWrapper(resource_stream(__name__, 'data/test_from_EcoliK12DH10B.fasta'))
        simulate_reads(inputfastafile,readlength=150,outfile=outputfasta)
        self.assertEqual(hashlib.sha224(outputfasta.getvalue().encode('latin-1')).hexdigest(),'91f0a1378bdf919d8e731237183950e3a8106231b11f7bd9832db600')
        pass
    
if __name__ == '__main__':
    unittest.main()