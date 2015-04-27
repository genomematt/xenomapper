#!/usr/bin/env python3
# encoding: utf-8
"""
test_mappability.py

Created by Matthew Wakefield on 2012-08-16.
Copyright (c) 2012 Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""

import unittest
import sys, io
from xenomapper.mappability import *
import hashlib
from pkg_resources import resource_stream
from string import ascii_uppercase, ascii_lowercase


__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2015 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "1.0b4"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Beta"

class test_mappability(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_parse_fasta(self):
        testdata = io.StringIO('>firstsequence\nGACAT\n>secondsequence\nGNATCAT')
        self.assertEqual(list(parse_fasta(testdata)),[('firstsequence', 'GACAT'), ('secondsequence', 'GNATCAT')])
        pass
    
    def test_simulate_reads(self):
        outputfasta = io.StringIO()
        inputfastafile = io.TextIOWrapper(resource_stream(__name__, 'data/test_from_EcoliK12DH10B.fasta'))
        simulate_reads(inputfastafile,readlength=150,outfile=outputfasta)
        self.assertEqual(hashlib.sha224(outputfasta.getvalue().encode('latin-1')).hexdigest(),'91f0a1378bdf919d8e731237183950e3a8106231b11f7bd9832db600')

        canned_result = '>testing_1\nabcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWX\n' + \
                        '>testing_2\nbcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXY\n' + \
                        '>testing_3\ncdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ\n'
        testdata = io.StringIO('>testing\n'+ascii_lowercase+ascii_uppercase)
        testout = io.StringIO('')
        simulate_reads(testdata,readlength=50,outfile=testout)
        self.assertEqual(canned_result, testout.getvalue())

        pass
    
    def test_make_blocklist(self):
        canned_result = ['abcdefghij', 'klmnopqrst', 'uvwxyzABCD',
                        'EFGHIJKLMN', 'OPQRSTUVWX', 'YZ']
        result = make_blocklist(ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
        pass
    
    def test_slice_string_in_blocks(self):
        canned_result = 'abcdefghij\nklmnopqrst\nuvwxyzABCD\nEFGHIJKLMN\nOPQRSTUVWX\nYZ\n'
        result = slice_string_in_blocks(ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
        pass
    
    def test_format_fasta(self):
        canned_result = '>fooGene\nabcdefghij\nklmnopqrst\nuvwxyzABCD\nEFGHIJKLMN\nOPQRSTUVWX\nYZ\n'
        result = format_fasta('fooGene',ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
        pass
    
    def test_smoothed_list(self):
        testdata = [1,2,3]*10 + [100,] + [1,2,3]*10
        result = smoothed_list(testdata)
        self.assertEqual(sum(result),219.8090909090909)
        self.assertEqual(max(result),6.714285714285714)
        self.assertEqual(min(result),1.9)
        self.assertEqual(len(result),len(testdata))
        pass
    
    def test_normalised_list(self):
        testdata = [1,2,3]*10 + [10,] + [1,2,3]*10
        result = normalised_list(testdata)
        self.assertAlmostEqual(sum(result),1.0)
        self.assertEqual(max(result),0.07692307692307693)
        self.assertEqual(min(result),0.007692307692307693)
        self.assertEqual(result[0]/result[2],1/3)
        self.assertEqual(len(result),len(testdata))
        pass
    
    def test_remove_small_values(self):
        testdata = list(range(100))
        result = remove_small_values(testdata)
        self.assertEqual(result[:10],[0,]*10)
        self.assertEqual(result[10:],testdata[10:])
        self.assertEqual(len(result),len(testdata))
        pass
    
    def test_mate_distribution_from_sam(self):
        hssam = io.TextIOWrapper(resource_stream(__name__, 'data/paired_end_testdata_human.sam'))
        result = mate_distribution_from_sam(samfile=hssam,sample_size=3)
        self.assertEqual(result, [0.0,]*164 + [0.047619047619047596,]*21 + [0.0,]*266)
        hssam.close()
        pass
    
    def test_single_end_mappability_from_sam(self):
        EcoliK12DH10B_150sam = io.TextIOWrapper(resource_stream(__name__, 'data/test_from_EcoliK12DH10B_150reads.sam'))
        resultfile = io.StringIO()
        single_end_mappability_from_sam(EcoliK12DH10B_150sam, outfile=resultfile, chromosome_sizes={'Chromosome':2752,'Repeat':991})
        self.assertEqual(hashlib.sha224(resultfile.getvalue().encode('latin-1')).hexdigest(),'af1a1a05670a6137aebb248e60173eb10275459318f3d0747bc5bf34')
        EcoliK12DH10B_150sam.close()
        resultfile.close()
        pass
    
    def test_mappability_obj(self):
        mappable = Mappability(chromosome_sizes={'Chromosome':20,})
        self.assertEqual(mappable['Chromosome'],[0,]*20)
        pair_mappability = mappable.single_end_to_paired(mate_density = [0,0.5,0.5,])
        self.assertEqual(pair_mappability['Chromosome'],[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0])
        mappable['Y'] = [0.0,]*15

        resultfile = io.StringIO()
        resultfile2 = io.StringIO()
        
        mappable.to_wiggle(wigglefile=resultfile, chromosomes = ['Chromosome'])
        self.assertEqual(resultfile.getvalue(),'fixedStep\tchrom=Chromosome\tstart=1\tstep=1\n' + '0\n'*20)
        resultfile.close()
        
        mappable = Mappability(chromosome_sizes={'X':30,})
        mappable['X'] = [0,1,1]*10
        mappable.to_wiggle(wigglefile=resultfile2)
        self.assertEqual(resultfile2.getvalue(),'fixedStep\tchrom=X\tstart=1\tstep=1\n' + '0\n1\n1\n'*10)
        resultfile2.close()
        
        pair_mappability = mappable.single_end_to_paired(mate_density = [0,0.4,0.5,0.1])
        self.assertEqual(pair_mappability['X'],[0.9, 1.0, 1.0,]*10 )
        pass
    
    def test_paired_end_mappability(self):
        resultfile = io.StringIO()
        mate_density = [0,0,0,0,0,0,0,0,0,0.01,0.45,0.41,0.13,0,0,0,0]
        wigglefile = io.StringIO('fixedStep\tchrom=Chromosome\tstart=1\tstep=1\n' + '1\n0\n'*50 + 'fixedStep\tchrom=Repeat\tstart=1\tstep=1\n' + '0\n'*10)
        paired_end_mappability(wigglefile,mate_density, outfile=resultfile,chromosome_sizes= {'Chromosome':100,'X':10})
        self.assertEqual(resultfile.getvalue(),'fixedStep\tchrom=Chromosome\tstart=1\tstep=1\n'+ '1.0\n0.42\n'*44 + '1.0\n0.01\n' + '1.0\n0.0\n'*5 + 'fixedStep\tchrom=X\tstart=1\tstep=1\n' + '0.0\n'*10)
        wigglefile.close()
        resultfile.close()
    
if __name__ == '__main__':
    unittest.main()