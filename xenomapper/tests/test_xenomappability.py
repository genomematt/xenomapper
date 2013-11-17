#!/usr/bin/env python3
# encoding: utf-8
"""
test_xenomappability.py

Created by Matthew Wakefield on 2012-08-16.
Copyright (c) 2012 Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""

import unittest
import sys, io
from xenomapper.xenomappability import simulate_reads
import hashlib
from pkg_resources import resource_stream

class test_main(unittest.TestCase):
    def setUp(self):
        pass
    def test_simulate_reads(self):
        outputfasta = io.StringIO()
        inputfastafile = io.TextIOWrapper(resource_stream(__name__, 'data/test_from_EcoliK12DH10B.fasta'))
        simulate_reads(inputfastafile,readlength=150,outfile=outputfasta)
        self.assertEqual(hashlib.sha224(outputfasta.getvalue().encode('latin-1')).hexdigest(),'3f9429eb9e50f90f7d5b79705d9e7ef5dabd90f1dcbd27c022650109')
        pass
    
if __name__ == '__main__':
    unittest.main()