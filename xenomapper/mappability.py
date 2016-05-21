#!/usr/bin/env python3
# encoding: utf-8
"""
mappability.py

This is an experimental addition to the Xenomapper tool for parsing pairs of sam files and returning 
sam files containing only reads where no better mapping exist in other files.
Used for filtering reads where multiple species may contribute (eg human tissue xenografted into mouse).

This file contains functionality for calculating paired end mappability estimates using a distribution
of insert sizes observed in the data and the probabilistic inference of a pair rescuing an otherwise
non-unique read.  This process is computationally intensive and in practice only improves mapping as
very low coverages.

This module should be considered experimental, and unlike the core xenomapper code which has been in
active use for several years, should be considered beta quality.  Caveat Emptor.

Created by Matthew Wakefield on 2011-12-08.
Copyright (c) 2011-2016  Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""
import sys
import os
import argparse
from statistics import *
from collections import Counter
from xenomapper.xenomapper import get_sam_header

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2016 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "1.0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development"


class Mappability(dict):
    def __init__(self, chromosome_sizes = {}):
        if chromosome_sizes:
            for chrom in chromosome_sizes:
                self[chrom] = [0,]*chromosome_sizes[chrom]
        self.chromosome_sizes = chromosome_sizes
        pass
    
    def to_wiggle(self, wigglefile=sys.stdout, chromosomes=[]):
        """Output mappability data to file in wiggle format"""
        ## Wiggle file format is:
        #fixedStep chrom=chrN start=pos step=1
        #value
        #value
        for chrom in sorted(self):
            if not chromosomes or chrom in chromosomes:
                print('fixedStep\tchrom={0}\tstart=1\tstep=1'.format(chrom), file=wigglefile)
                for score in self[chrom]:
                    print(str(score), file=wigglefile)
        pass
    
    def from_wiggle(self,wigglefile=sys.stdin,datatype=float):
        """Load mapability data from a wiggle file
        The wiggle file must be fixed step format, step of one and start at 1
        Multiple chromosomes may be present in the same file
        This function will overwrite any existing chromosomes
        with the same name but can be used sequentially for
        different chromosomes
        """
        chrom=None
        values=[]
        for line in wigglefile:
            if line.startswith('fixedStep'):
                #add a chromosome to self
                #check it is not already there
                #check start=1 and step=1
                line = line.strip().split('\t')
                if line[1].split('=')[0] != 'chrom' or line[2] != 'start=1' or line[3] != 'step=1': #pragma: no cover
                    raise ValueError('Unsupported wiggle fixed step format [must be in the format "fixedStep chrom=chrX start=1 step=1"] {0}'.format(line))
                if chrom and values:
                    #if we have accumulated values for a previous chromosome add them to self
                    self[chrom]=values
                chrom=line[1].split('=')[1]
                values=[]
            else:
                try:
                    values.append(datatype(line))
                except: #pragma: no cover
                    #placeholder for type conversion error handling
                    raise
        #at end of file add remaining data to self
        self[chrom]=values
        for chrom in self:
            self.chromosome_sizes[chrom]=len(self[chrom])
        pass
    
    def single_end_to_paired(self, mate_density = [1,]):
        """Produce a new mappability object with paired end mapping probilities
        Defines paired end mappability as either end being uniquely mappable.
        Arguments:
            mate_density:   a list of floats between 0.0 and 1.0 representing mate densities.
                            First entry corresponds to current position
                            all entries must sum to 1.0
        """
        def _mappability_by_mate_density(mapability,mate_density):
            #defined as local scope function as may be replaced for speed.
            result = 0.0
            j = 0
            while j < len(mapability) and j < len(mate_density):
                result += mapability[j] * mate_density[j]
                j += 1
            return result
        
        #should sum to 1 but allow for numerical error
        #summing to 1 is not algorythmically essential
        assert abs(sum(mate_density)-1.0) < 0.000001 
        
        paired_mappability = Mappability(chromosome_sizes = self.chromosome_sizes)
        
        for chrom in self:
            for i in range(len(self[chrom])):
                if self[chrom][i] == 1:
                    paired_mappability[chrom][i] = 1.0
                else:
                    paired_mappability[chrom][i] = _mappability_by_mate_density(self[chrom][i:i+len(mate_density)],mate_density)
        
        return paired_mappability
    

def parse_fasta(fastafile, token='>'):
    """fasta and multi-fasta file parser
    Usage: for name,seq in fasta(open(filename)):
                do something
           next(parse_fasta(open(filename)))
    """
    with fastafile as f:
        seq = None
        name = None   
        for line in f:
            line = line.strip()
            if line.startswith(token):
                if name:
                    yield (name, seq)
                seq = ''
                name = line[1:]
            elif seq != None:
                seq += line
        if name:
            yield (name, seq)

def make_blocklist(seqstring, block_size=80):
    """format sequence into a list of blocks"""
    blocklist = []
    seqlength = len(seqstring)
    for block in range(0, seqlength, block_size):
        if block + block_size < seqlength:
            blocklist.append(seqstring[block: block + block_size])
        else:
            blocklist.append(seqstring[block:])
    return blocklist

def slice_string_in_blocks(seqstring, block_size=80):
    """slice string into block_size[=80] lines"""
    blocklist = make_blocklist(seqstring, block_size)
    return '\n'.join(blocklist) + '\n'

def format_fasta(name,seq, block_size=80):
    """returns a string in fasta format"""
    return '>'+name+'\n'+slice_string_in_blocks(seq,block_size)

def simulate_reads(fastafile, readlength=100, outfile=sys.stdout):
    for name, seq in parse_fasta(fastafile):
        for x in range(len(seq)-readlength+1):
            newname = '{chrom}_{one_based_pos}'.format(chrom=name.split()[0],one_based_pos=x+1)
            outfile.write(format_fasta(newname,seq[x:x+readlength]))
    pass

def single_end_mappability_from_sam(samfile, outfile=sys.stdout, fill_sequence_gaps=True, chromosome_sizes = {}):
    mappable = Mappability(chromosome_sizes=chromosome_sizes)
    
    #parse sam data and add to object
    header = get_sam_header(samfile)
    mappable_chrom = None
    mappable_pos = 0
    mappable_values = []
    for line in samfile:
        name, x, chrom, pos, flag, cigar, mate_chr, mate_pos, insert_size, seq, qual, *tags = line.strip('\n').split()
        true_chrom = '_'.join(name.split('_')[:-1]) #name may have _ so we split off last and join
        if mappable_chrom != true_chrom: #add next chromosome to mappable_chrom if new
            if mappable_chrom and mappable_values:
                mappable[mappable_chrom] = mappable_values
            mappable_chrom = chrom
            mappable_pos = 0
            mappable_values = []
        mappable_pos += 1
        name_pos = int(name.split('_')[-1])
        if name_pos != mappable_pos: #pragma: no cover
            if not name_pos > mappable_pos:
                raise ValueError('Name is not sequential.  SAM must be in name sorted order Name: {0} Expected: {1}_{2}'.format(name,mappable_chrom,mappable_pos))
            else:
                #there are missing reads in the sequence.
                number_missing = name_pos - mappable_pos
                mappable_pos = name_pos
                mappable_values.extend([0,]*number_missing)
                
        if (true_chrom, name_pos) == (chrom, int(pos)) and flag == '42':
            mappable_values.append(1)
        else:
            mappable_values.append(0)
    #add remainging chromosome and values at the end of the file
    if mappable_chrom and mappable_values:
        mappable[mappable_chrom] = mappable_values
    
    #write wiggle file
    mappable.to_wiggle(wigglefile=outfile)
    
    pass

def paired_end_mappability(wiggle, mate_density, outfile=sys.stdout, chromosome_sizes={}):
    """Create a wiggle file of read mappability using inferred mapping rate of pairs
    Arguments:
        wiggle           - a wiggle file of mappability
        mate_density     - an iterable of floats summing to 1 representing
                           the probability of observing a pair.
                           (Usually output of mate_distribution_from_sam)
        outfile          - file object for writing ouput wiggle file
        chromosome_sizes - a dictionary of chromosome sizes
    """
    mappable = Mappability(chromosome_sizes=chromosome_sizes)
    chromosomes = list(chromosome_sizes.keys())
    
    mappable.from_wiggle(wiggle, datatype=float)
    
    pair_mappability = mappable.single_end_to_paired(mate_density = mate_density)
    
    pair_mappability.to_wiggle(wigglefile=outfile, chromosomes=chromosomes)
    
    pass

def smoothed_list(the_list,width=10):
    return [mean(the_list[max(0,x-width-1):x+width]) for x in range(len(the_list))]

def normalised_list(the_list):
    total = sum(the_list)
    return [x/total for x in the_list]

def remove_small_values(the_list,relative_limit=0.1):
    """replaces values less than maximum_value * relative_limit with 0"""
    min_value = max(the_list) * relative_limit
    result = []
    for x in the_list:
        if x > min_value:
            result.append(x)
        else:
            result.append(0)
    return result
    
def mate_distribution_from_sam(samfile=sys.stdin, sample_size=10000):
    """Calculate the mate density (distribution) of read pair sizes from a sam file"""
    sizes = []
    for line in samfile:
        if not line or line[0] == '@':
            continue
        name, x, chrom, pos, flag, cigar, mate_chr, mate_pos, insert_size, seq, qual, *tags = line.strip('\n').split()
        if flag not in [] and insert_size != '0':#this should be for flag in and list of valid mate pair mapped flags
            sizes.append(abs(int(insert_size)))
        if sample_size and len(sizes) > sample_size:
            break
    frequencies = Counter(sizes)
    mate_density = []
    for i in range(0,max(sizes)):
        if i in frequencies:
            mate_density.append(frequencies[i])
        else:
            mate_density.append(0)
    return normalised_list(remove_small_values(smoothed_list(mate_density)))

def command_line_interface(): #pragma: no cover
    parser = argparse.ArgumentParser(description='Caution: Experimental - Beta quality functionality.\
                                                A script for generating mappability estimates for paired end data.\
                                                Paired end mappability is inferred from single end mappability.\
                                                Step one is to generate a fasta file of reads from a fasta file using --fasta \
                                                You then process these reads through your mapping process. \
                                                For xenomapper this is mapping against the two reference genomes in \
                                                single end mode, and processing to primary unique mappings. \
                                                Step three is to generate a single end mappability wiggle with --mapped_test_data.\
                                                Step four is to generate a mappability wiggle file using the --single_end_wiggle and --sam_for_sizes.')
    parser.add_argument('--fasta',
                        type=argparse.FileType('rt'),
                        help='Process a fasta genome file of sequences to simulated reads in fasta format.\
                                Outputs fasta file to standard output.')
    parser.add_argument('--readlength',
                        type=int,
                        default = 100,
                        help='The readlength to simulate.')
    parser.add_argument('--mapped_test_data',
                        type=argparse.FileType('rt'),
                        help='a SAM format input file of mapped reads generated by the --fasta command.\
                              The SAM file must be sorted by read name.\
                              outputs a fixed step wiggle file to standard output.')
    parser.add_argument('--single_end_wiggle',
                        type=argparse.FileType('rt'),
                        help='a wiggle file of single end mappabilities')
    parser.add_argument('--sam_for_sizes',
                        type=argparse.FileType('rt'),
                        help='a sam file for calculating insert sizes')
    parser.add_argument('--version',
                        action='store_true',
                        help='print version information and exit')
    args = parser.parse_args()
    if args.version:
        print(__version__)
        sys.exit()
    if (not args.fasta) and (not args.mapped_test_data) and (not args.single_end_wiggle):
        print('ERROR: Insufficient arguments provided')
        parser.print_help()
        sys.exit(1)
    return args

def main(args=None): #pragma: no cover
    if not args:
        args = command_line_interface()
    if args.fasta:
        simulate_reads(fastafile=args.fasta, readlength=args.readlength)
    elif args.mapped_test_data:
        single_end_mappability_from_sam(samfile=args.mapped_test_data)
    elif args.single_end_wiggle:
        if not args.sam_for_sizes:
            raise RuntimeError('You must provide a sam file to estimate the mate pair distance distribution')
        mate_density = mate_distribution_from_sam(args.sam_for_sizes)
        paired_end_mappability(wiggle=args.single_end_wiggle, mate_density=mate_density)
    pass


if __name__ == '__main__': # pragma: no cover
    main()
    
