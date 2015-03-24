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