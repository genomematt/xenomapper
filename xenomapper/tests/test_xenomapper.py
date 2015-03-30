#!/usr/bin/env python3
# encoding: utf-8
"""
test_xenomapper.py

Created by Matthew Wakefield on 2012-08-16.
Copyright (c) 2012 Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""

import unittest
import sys, io
from xenomapper.xenomapper import *
import hashlib
from pkg_resources import resource_stream

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2015 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development"

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
        
    def test_get_mapping_state(self):
        inpt_and_outpt = [
                            ((200,199,199,198,float('-inf')),'primary_specific'),
                            ((200,200,199,198,float('-inf')),'primary_multi'),
                            ((199,198,200,198,float('-inf')),'secondary_specific'),
                            ((199,198,200,200,float('-inf')),'secondary_multi'),
                            #((float('-inf'),float('-inf'),float('-inf'),float('-inf'),float('-inf')),'unassigned'), #need to test vs minimum & change to <=
                            ((200,199,200,198,float('-inf')),'unresolved'),
                            ((200,199,199,199,float('-inf')),'primary_specific'),
                            ((200,200,199,199,float('-inf')),'primary_multi'),
                            ((199,199,200,199,float('-inf')),'secondary_specific'),
                            ((199,199,200,200,float('-inf')),'secondary_multi'),
                            ((9,8,8,8,10),'unassigned'),
                            ((200,200,200,200,float('-inf')),'unresolved'),
                            ((-6,float('-inf'),float('-inf'),float('-inf'),float('-inf')),'primary_specific'),
                            ((float('-inf'),float('-inf'),-6,float('-inf'),float('-inf')),'secondary_specific'),
                            ((-6,float('-inf'),-2,float('-inf'),float('-inf')),'secondary_specific'),
                            ((0,float('-inf'),-2,float('-inf'),float('-inf')),'primary_specific'),
                            ((-2,float('-inf'),0,float('-inf'),float('-inf')),'secondary_specific'),
                           
        ]
        
        for inpt, outpt in inpt_and_outpt:
            self.assertEqual(get_mapping_state(*inpt),outpt)
        pass
    
    def test_get_cigarbased_AS_tag(self):
        inpt_and_outpt = [
                        (['HWI-ST960:63:D0CYJACXX:4:1101:21264:2228', '4', '*', '0', '0', '*', '*', '0', '0', 'TGGTAGTATTGGTTATGGTTCATTGTCCGGAGAGTATATTGTTGAAGAGG', 'BBCBDFDDHHHGFHHIIIIIJIJJJIGJJJGIAF:CFEGHGGHEEEG@HI', 'YT:Z:UU'],
                        float('-inf')),
                        (['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0'],0),
                        (['', '', '', '', '', '1S49M', '', '', '', '', '', 'NM:i:0'],-2),
                        (['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:2'],-12),
                        (['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:i:99'],0),
                        (['', '', '', '', '', '10M1I39M', '', '', '', '', '', 'NM:i:0'],-8),
                        (['', '', '', '', '', '10M1D39M', '', '', '', '', '', 'NM:i:0'],-8),
                        (['', '', '', '', '', '10M2D38M', '', '', '', '', '', 'NM:i:0'],-11),
                        (['', '', '', '', '', '10M1I10M1D28M', '', '', '', '', '', 'NM:i:0'],-16),
                        (['', '', '', '', '', '10M1234N40M', '', '', '', '', '', 'NM:i:0'],0),
                        ]
        
        for inpt, outpt in inpt_and_outpt:
            self.assertEqual(get_cigarbased_AS_tag(inpt),outpt)
        pass

    
if __name__ == '__main__':
    unittest.main()