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

    def test_parse_fasta(self):
        testdata = io.StringIO('>firstsequence\nGACAT\n>secondsequence\nGNATCAT')
        self.assertEqual(list(parse_fasta(testdata)),[('firstsequence', 'GACAT'), ('secondsequence', 'GNATCAT')])
    
    def test_make_blocklist(self):
        canned_result = ['abcdefghij', 'klmnopqrst', 'uvwxyzABCD',
                        'EFGHIJKLMN', 'OPQRSTUVWX', 'YZ']
        result = make_blocklist(ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
    
    def test_slice_string_in_blocks(self):
        canned_result = 'abcdefghij\nklmnopqrst\nuvwxyzABCD\nEFGHIJKLMN\nOPQRSTUVWX\nYZ\n'
        result = slice_string_in_blocks(ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
    
    def test_format_fasta(self):
        canned_result = '>fooGene\nabcdefghij\nklmnopqrst\nuvwxyzABCD\nEFGHIJKLMN\nOPQRSTUVWX\nYZ\n'
        result = format_fasta('fooGene',ascii_lowercase+ascii_uppercase, 10)
        self.assertEqual(canned_result, result)
    
    def test_simulate_reads(self):
        canned_result = '>testing_1\nabcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWX\n' + \
                        '>testing_2\nbcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXY\n' + \
                        '>testing_3\ncdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ\n'
        testdata = io.StringIO('>testing\n'+ascii_lowercase+ascii_uppercase)
        testout = io.StringIO('')
        simulate_reads(testdata,readlength=50,outfile=testout)
        self.assertEqual(canned_result, testout.getvalue())
    
    def test_smoothed_list(self):
        testdata = [1,2,3]*10 + [100,] + [1,2,3]*10
        result = smoothed_list(testdata)
        self.assertEqual(sum(result),219.8090909090909)
        self.assertEqual(max(result),6.714285714285714)
        self.assertEqual(min(result),1.9)
        self.assertEqual(len(result),len(testdata))
    
    def test_normalised_list(self):
        testdata = [1,2,3]*10 + [10,] + [1,2,3]*10
        result = normalised_list(testdata)
        self.assertAlmostEqual(sum(result),1.0)
        self.assertEqual(max(result),0.07692307692307693)
        self.assertEqual(min(result),0.007692307692307693)
        self.assertEqual(result[0]/result[2],1/3)
        self.assertEqual(len(result),len(testdata))
    
    def test_remove_small_values(self):
        testdata = list(range(100))
        result = remove_small_values(testdata)
        self.assertEqual(result[:10],[0,]*10)
        self.assertEqual(result[10:],testdata[10:])
        self.assertEqual(len(result),len(testdata))

#def mate_distribution_from_sam(samfile=sys.stdin, sample_size=10000):
#    sizes = []
#    for line in samfile:
#        if not line or line[0] == '@':
#            continue
#        name, x, chrom, pos, flag, cigar, mate_chr, mate_pos, insert_size, seq, qual, *tags = line.strip('\n').split()
#        if flag not in [] and insert_size != '0':#this should be for flag in and list of valid mate pair mapped flags
#            sizes.append(abs(int(insert_size)))
#        if sample_size and len(sizes) > sample_size:
#            break
#    frequencies = Counter(sizes)
#    mate_density = []
#    for i in range(0,max(sizes)):
#        if i in frequencies:
#            mate_density.append(frequencies[i])
#        else:
#            mate_density.append(0)
#    return normalised_list(remove_small_values(smoothed_list(mate_density)))

#def single_end_mappability_from_sam(samfile, outfile=sys.stdout, fill_sequence_gaps=True, chromosome_sizes = {}):
#    mapable = Mappability(chromosome_sizes=chromosome_sizes)
#    
#    #parse sam data and add to object
#    header = get_sam_header(samfile)
#    mapable_chrom = None
#    mapable_pos = 0
#    mapable_values = []
#    for line in samfile:
#        name, x, chrom, pos, flag, cigar, mate_chr, mate_pos, insert_size, seq, qual, *tags = line.strip('\n').split()
#        if mapable_chrom != name.split('_')[0]:
#            if mapable_chrom and mapable_values:
#                mapable[mapable_chrom] = mapable_values
#            mapable_chrom = chrom
#            mapable_pos = 0
#            mapable_values = []
#        mapable_pos += 1
#        name_pos = int(name.split('_')[1])
#        if name_pos != mapable_pos:
#            if not name_pos > mapable_pos:
#                raise ValueError('Name is not sequential.  SAM must be in name sorted order Name: {0} Expected: {1}_{2}'.format(name,mapable_chrom,mapable_pos))
#            else:
#                #there are missing reads in the sequence.
#                number_missing = name_pos - mapable_pos
#                mapable_pos = name_pos
#                mapable_values.extend([0,]*number_missing)
#                
#        if name.split('_') == [chrom, pos] and flag == '42':
#            mapable_values.append(1)
#        else:
#            mapable_values.append(0)
#    #add remainging chromosome and values at the end of the file
#    if mapable_chrom and mapable_values:
#        mapable[mapable_chrom] = mapable_values
#    
#    #write wiggle file
#    mapable.to_wiggle(wigglefile=outfile)
#    
#    pass
#
#def paired_end_mappability(wiggle, mate_density, outfile=sys.stdout, chromosome_sizes={}):
#    mapable = Mappability(chromosome_sizes=chromosome_sizes)
#    
#    mapable.from_wiggle(wiggle, datatype=int)
#    
#    pair_mapability = mapable.single_end_to_paired(mate_density = mate_density)
#    
#    pair_mapability.to_wiggle(wigglefile=outfile)
#    
#    pass
    
    
if __name__ == '__main__':
    unittest.main()