#!/usr/bin/env python3
# encoding: utf-8
"""
test_xenomapper.py

Created by Matthew Wakefield on 2012-08-16.
Copyright (c) 2011-2016 Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""

import unittest
import sys, io
from xenomapper.xenomapper import *
import hashlib
from pkg_resources import resource_stream

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2016 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production/Stable"

class test_main(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_process_headers(self):
        test_primary_specific_outfile = io.StringIO()
        test_secondary_specific_outfile = io.StringIO()
        test_primary_multi_outfile = io.StringIO()
        test_secondary_multi_outfile = io.StringIO()
        test_unassigned_outfile = io.StringIO()
        test_unresolved_outfile = io.StringIO()
        sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
        sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_mouse.sam'))
        process_headers(sam1,sam2,
                     primary_specific=test_primary_specific_outfile,
                     secondary_specific=test_secondary_specific_outfile,
                     primary_multi=test_primary_multi_outfile,
                     secondary_multi=test_secondary_multi_outfile,
                     unresolved=test_unresolved_outfile,
                     unassigned=test_unassigned_outfile,
                     )
        self.assertEqual(len(test_primary_specific_outfile.getvalue()),695)
        self.assertEqual(len(test_secondary_specific_outfile.getvalue()),629)
        self.assertEqual(len(test_primary_multi_outfile.getvalue()),708)
        self.assertEqual(len(test_secondary_multi_outfile.getvalue()),642)
        self.assertEqual(len(test_unassigned_outfile.getvalue()),705)
        self.assertEqual(len(test_unresolved_outfile.getvalue()),705)
        sam1.close()
        sam2.close()
        pass

    def test_consistent_output_SE(self):
        test_primary_specific_outfile = io.StringIO()
        test_secondary_specific_outfile = io.StringIO()
        test_primary_multi_outfile = io.StringIO()
        test_secondary_multi_outfile = io.StringIO()
        test_unassigned_outfile = io.StringIO()
        test_unresolved_outfile = io.StringIO()
        sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/test_human_in.sam'))
        sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/test_mouse_in.sam'))
        process_headers(sam1,sam2,
                         primary_specific=test_primary_specific_outfile,
                         secondary_specific=test_secondary_specific_outfile,
                         primary_multi=test_primary_multi_outfile,
                         secondary_multi=test_secondary_multi_outfile,
                         unresolved=test_unresolved_outfile,
                         unassigned=test_unassigned_outfile,
                         )
        cat_counts = main_single_end(getReadPairs(sam1,sam2),
                         primary_specific=test_primary_specific_outfile,
                         secondary_specific=test_secondary_specific_outfile,
                         primary_multi=test_primary_multi_outfile,
                         secondary_multi=test_secondary_multi_outfile,
                         unresolved=test_unresolved_outfile,
                         unassigned=test_unassigned_outfile,
                         )
        self.assertEqual(cat_counts['primary_specific'],
                         len(test_primary_specific_outfile.getvalue().split('\n'))-4) #29 lines of header in this file
        self.assertEqual(cat_counts['primary_multi'],
                         len(test_primary_multi_outfile.getvalue().split('\n'))-4) #29 lines of header in this file
        self.assertEqual(cat_counts['secondary_specific'],
                         len(test_secondary_specific_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
        self.assertEqual(cat_counts['secondary_multi'],
                         len(test_secondary_multi_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
        self.assertEqual(cat_counts['unassigned'],
                         len(test_unassigned_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
        self.assertEqual(cat_counts['unassigned'],
                         len(test_unassigned_outfile.getvalue().split('\n'))-4) #26 lines of header in this file
        self.assertEqual(hashlib.sha224(test_primary_specific_outfile.getvalue().encode('latin-1')).hexdigest(),'3512efc4f7b0e37bbcd2871826628ea138112086880461bc2c40cc93')
        sam1.close()
        sam2.close()
        pass

    def test_consistent_output_PE(self):
        test_primary_specific_outfile = io.StringIO()
        test_secondary_specific_outfile = io.StringIO()
        test_primary_multi_outfile = io.StringIO()
        test_secondary_multi_outfile = io.StringIO()
        test_unassigned_outfile = io.StringIO()
        test_unresolved_outfile = io.StringIO()
        sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
        sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_mouse.sam'))
        process_headers(sam1,sam2,primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
        cat_counts = main_paired_end(getReadPairs(sam1,sam2),
                                     primary_specific=test_primary_specific_outfile,
                                     secondary_specific=test_secondary_specific_outfile,
                                     primary_multi=test_primary_multi_outfile,
                                     secondary_multi=test_secondary_multi_outfile,
                                     unresolved=test_unresolved_outfile,
                                     unassigned=test_unassigned_outfile,
                                     )
        self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'primary_specific' in x])*2,
                         len(test_primary_specific_outfile.getvalue().split('\n'))-30) #29 lines of header in this file
        self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'secondary_specific' in x and not 'primary_specific' in x])*2,
                         len(test_secondary_specific_outfile.getvalue().split('\n'))-27) #26 lines of header in this file
        self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'primary_multi' in x and not 'primary_specific' in x and not 'secondary_specific' in x])*2,
                        len(test_primary_multi_outfile.getvalue().split('\n'))-1)
        self.assertEqual(sum([cat_counts[x] for x in cat_counts if 'secondary_multi' in x \
                        and not 'primary_multi' in x and not 'primary_specific' in x and not 'secondary_specific' in x])*2,
                        len(test_secondary_multi_outfile.getvalue().split('\n'))-1)
        self.assertEqual(hashlib.sha224(test_primary_specific_outfile.getvalue().encode('latin-1')).hexdigest(),'ecba2de3e3af9c7405a84ad2a4ebaf194ebfb4df76f45c311c0f681d')
        sam1.close()
        sam2.close()
        pass

    def test_consistent_output_conservative_PE(self):
        test_primary_specific_outfile = io.StringIO()
        test_secondary_specific_outfile = io.StringIO()
        test_primary_multi_outfile = io.StringIO()
        test_secondary_multi_outfile = io.StringIO()
        test_unassigned_outfile = io.StringIO()
        test_unresolved_outfile = io.StringIO()
        sam1 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
        sam2 = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_mouse.sam'))
        process_headers(sam1,sam2,primary_specific=test_primary_specific_outfile, secondary_specific=test_secondary_specific_outfile)
        cat_counts = conservative_main_paired_end(getReadPairs(sam1,sam2),
                                                 primary_specific=test_primary_specific_outfile,
                                                 secondary_specific=test_secondary_specific_outfile,
                                                 primary_multi=test_primary_multi_outfile,
                                                 secondary_multi=test_secondary_multi_outfile,
                                                 unresolved=test_unresolved_outfile,
                                                 unassigned=test_unassigned_outfile,
                                                 )
        self.assertEqual(cat_counts[('primary_specific', 'secondary_specific')]*4 +\
                         cat_counts[('unresolved', 'unresolved')]*4, len(test_unresolved_outfile.getvalue().split('\n'))-1) #this test is not exhastive. Only states in test data
        self.assertEqual(cat_counts[('primary_multi', 'primary_multi')]*2, len(test_primary_multi_outfile.getvalue().split('\n'))-1)
        self.assertEqual(cat_counts[('secondary_multi', 'secondary_multi')]*2, len(test_secondary_multi_outfile.getvalue().split('\n'))-1)
        self.assertEqual(cat_counts[('primary_specific', 'primary_specific')]*2 + cat_counts[('primary_specific', 'primary_multi')]*2 + cat_counts[('primary_multi', 'primary_specific')]*2,
                         len(test_primary_specific_outfile.getvalue().split('\n'))-30) #29 lines of header in this file
        self.assertEqual(cat_counts[('secondary_specific', 'secondary_specific')]*2 + cat_counts[('secondary_specific', 'secondary_multi')]*2 + cat_counts[('secondary_multi', 'secondary_specific')]*2,
                         len(test_secondary_specific_outfile.getvalue().split('\n'))-27)  #26 lines of header in this file
        self.assertEqual(cat_counts[('unassigned', 'unassigned')]*2, len(test_unassigned_outfile.getvalue().split('\n'))-1)
        
        self.assertEqual(hashlib.sha224(test_primary_specific_outfile.getvalue().encode('latin-1')).hexdigest(),'8f5349ac96f194a4600bf0542cb1a6ebf71ada14b8ee0986598d7f58')
        sam1.close()
        sam2.close()
        pass

        
    def test_get_mapping_state(self):
        inpt_and_outpt = [
                            ((200,199,199,198,float('-inf')),'primary_specific'),
                            ((200,200,199,198,float('-inf')),'primary_multi'),
                            ((199,198,200,198,float('-inf')),'secondary_specific'),
                            ((199,198,200,200,float('-inf')),'secondary_multi'),
                            ((float('-inf'),float('-inf'),float('-inf'),float('-inf'),float('-inf')),'unassigned'),
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
    
    def test_get_tag(self):
        inpt_and_outpt = [
                        ((['HWI-ST960:63:D0CYJACXX:4:1101:21264:2228', '4', '*', '0', '0', '*', '*', '0', '0', 'TGGTAGTATTGGTTATGGTTCATTGTCCGGAGAGTATATTGTTGAAGAGG', 'BBCBDFDDHHHGFHHIIIIIJIJJJIGJJJGIAF:CFEGHGGHEEEG@HI', 'YT:Z:UU'],'AS'),
                        float('-inf')),
                        ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:101', 'XS:i:99'],'AS'),101),
                        ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:i:99'],'XS'),99),
                        ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:i:99'],'NM'),0),
                        ]
        for inpt, outpt in inpt_and_outpt:
            self.assertEqual(get_tag(*inpt),outpt)
        pass
    
    def test_get_tag_with_ZS_as_XS(self):
        inpt_and_outpt = [
                        ((['HWI-ST960:63:D0CYJACXX:4:1101:21264:2228', '4', '*', '0', '0', '*', '*', '0', '0', 'TGGTAGTATTGGTTATGGTTCATTGTCCGGAGAGTATATTGTTGAAGAGG', 'BBCBDFDDHHHGFHHIIIIIJIJJJIGJJJGIAF:CFEGHGGHEEEG@HI', 'YT:Z:UU'],'AS'),
                        float('-inf')),
                        ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:101', 'XS:A:+', 'ZS:i:99'],'AS'),101),
                        ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:A:+', 'ZS:i:99'],'XS'),99),
                        ((['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:A:+', 'ZS:i:99'],'NM'),0),
                        ]
        for inpt, outpt in inpt_and_outpt:
            self.assertEqual(get_tag_with_ZS_as_XS(*inpt),outpt)
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
        
        self.assertEqual(get_cigarbased_AS_tag(['', '', '', '', '', '50M', '', '', '', '', '', 'NM:i:0', 'AS:i:100', 'XS:i:99'],tag='XS'),99)
        pass
        
    def test_output_summary(self):
        test_outfile = io.StringIO()
        canned_output='--------------------------------------------------------------------------------\n' + \
                    'Read Count Category Summary\n\n' + \
                    '|       Category                                     |     Count       |\n' + \
                    '|:--------------------------------------------------:|:---------------:|\n' + \
                    '|  bar                                               |            101  |\n' + \
                    '|  foo                                               |              1  |\n\n'
        output_summary({'foo':1,'bar':101},outfile=test_outfile)
        self.assertEqual(test_outfile.getvalue(),canned_output)
        pass

    
if __name__ == '__main__':
    unittest.main()