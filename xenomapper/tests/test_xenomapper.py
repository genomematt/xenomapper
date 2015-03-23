#!/usr/bin/env python3
# encoding: utf-8
"""
test_xenomapper.py

Created by Matthew Wakefield on 2012-08-16.
Copyright (c) 2012 Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""

import unittest
import sys, io
from xenomapper.xenomapper import process_headers, getReadPairs, main_single_end, main_paired_end
import hashlib
from pkg_resources import resource_stream

class test_main(unittest.TestCase):
    def setUp(self):
        pass
    def test_consistent_output_SE(self):
        test_primary_specific_outfile = io.StringIO()
        test_secondary_specific_outfile = io.StringIO()
        test_primary_multi_outfile = io.StringIO()
        test_secondary_multi_outfile = io.StringIO()
        test_unassigned_outfile = io.StringIO()
        #sam1 = open('test_human_in.sam','r')
        #sam2 = open('test_mouse_in.sam','r')
        sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/test_human_in.sam'))
        sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/test_human_in.sam'))
        process_headers(sam1,sam2,primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
        main_single_end(getReadPairs(sam1,sam2), primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
        test_primary_specific_outfile.seek(0)
        self.assertEqual(hashlib.sha224(test_primary_specific_outfile.read().encode('latin-1')).hexdigest(),'f251316f5737d0bc4496ef4e695b14f553daef379dd5e7c6b456fdcf')
        sam1.close()
        sam2.close()
        pass

    def test_consistent_output_PE(self):
        test_primary_specific_outfile = io.StringIO()
        test_secondary_specific_outfile = io.StringIO()
        test_primary_multi_outfile = io.StringIO()
        test_secondary_multi_outfile = io.StringIO()
        test_unassigned_outfile = io.StringIO()
        #sam1 = open('paired_end_testdata_human.sam','r')
        #sam2 = open('paired_end_testdata_mouse.sam','r')
        sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
        sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_mouse.sam'))
        process_headers(sam1,sam2,primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
        main_paired_end(getReadPairs(sam1,sam2), primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
        test_primary_specific_outfile.seek(0)
        self.assertEqual(hashlib.sha224(test_primary_specific_outfile.read().encode('latin-1')).hexdigest(),'7b34b93efc8d8eb284fd294534f2ac82a85d41cf39039e477a227f2f')
        sam1.close()
        sam2.close()
        pass
    

    
if __name__ == '__main__':
    unittest.main()