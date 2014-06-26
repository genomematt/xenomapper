#!/usr/bin/env python3
# encoding: utf-8
"""
xenomapper.py

A script for parsing pairs of sam files and returning sam files containing only reads where no better mapping exist in other files.
Used for filtering reads where multiple species may contribute (eg human tissue xenografted into mouse).

Created by Matthew Wakefield on 2011-12-08.
Copyright (c) 2011  Matthew Wakefield and The Walter and Eliza Hall Institute. All rights reserved.
"""
import sys
import os
import argparse

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011,  Matthew Wakefield and The Walter and Eliza Hall Institute"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.1"
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
        for chrom in self:
            if chromosomes and not chrom in chromosomes:
                continue
            print('fixedStep\tchrom={0}\tstart=1\tstep=1'.format(chrom), file=wigglefile)
            for score in self[chrom]:
                print(str(score))
        pass
    
    def from_wiggle(self,wigglefile=sys.stdin):
        """Load mapability data from a wiggle file"""
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
            result = 0
            j = 0
            while j < len(mapability) and j < len(mate_density):
                result += mapability[j] * mate_density[j]
                j += 1
            return result
        
        assert sum(mate_density) == 1.0
        
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
           parse_fasta(open(filename)).next()
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
        
def command_line_interface():
    parser = argparse.ArgumentParser(description='A script for generating mappability estimates for paired end data.\
                                                Paired end mappability is inferred from single end mappability.\
                                                Step one is to generate a fasta file of reads from a fasta file using --fasta \
                                                You then process these reads through your mapping process. \
                                                For xenomapper this is mapping against the two reference genomes in \
                                                single end mode, and processing to primary unique mappings. \
                                                Step three is to generate a single end mappability vector with --mapped_test_data.\
                                                Step four is to generate a mappability bed file using the --binary_mapability and --bam_for_size_distribution.')
    parser.add_argument('--fasta',
                        type=argparse.FileType('rt'),
                        help='Process a fasta genome file of sequences to simulated reads in fasta format.\
                                Outputs fasta file to standard output.')
    parser.add_argument('--readlength',
                        type=int,
                        default = 100,
                        help='The readlength to simulate.')
    parser.add_argument('--mapped_test_data',
                        type=argparse.FileType('wt'),
                        help='a SAM format input file of mapped reads generated by the --fasta command.\
                              outputs a fixed step wiggle file to standard output.')
    return parser.parse_args()

def main():
    args = command_line_interface()
    if args.fasta:
        simulate_reads(fastafile=args.fasta, readlength=args.readlength)
    pass


if __name__ == '__main__':
    main()
    
