#!/usr/bin/env python3
# encoding: utf-8
"""
xenomapper.py

A script for parsing pairs of sam files and returning sam files containing only reads where no better mapping exist in other files.
Used for filtering reads where multiple species may contribute (eg human tissue xenografted into mouse).

Created by Matthew Wakefield.
Copyright (c) 2011-2014  Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne. All rights reserved.

   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


"""
import sys
import os
import argparse, textwrap
import subprocess

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2014 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.3.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development"

def get_sam_header(samfile):
    line = "@"
    header = []
    pointer = 0
    while line[0] == '@':
        pointer = samfile.tell()
        line = samfile.readline().strip('\n')
        if line[0] == '@':
            header.append(line)
    samfile.seek(pointer) #set file to first line after header
    return header

def get_bam_header(bamfile):
    p = subprocess.Popen('samtools view -H -',stdin=bamfile,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    header = []
    for line in p.stdout:
        header.append(line.decode('ascii'))
    bamfile.seek(0) #reset to start of file for next samtools call
    return [x.strip('\n') for x in header]

    
#### To do: Check aligner and command line options in SAM file & warn if not bowtie2 or if there is an argument mismatch


#['HWI-ST960:63:D0CYJACXX:4:1101:6951:2219', '0', '9', '20953017', '44', '49M1S', '*', '0', '0', 'GNTTTATTGGAGCAGCTATTGGCTTCTTCATTGCAGGAGGAAAAAAAGGT', '@#1ADDDFFHFFHIJIJJIJIJJJJ?EH@@FFFGIII@GHCHIIJJI@F8', 'AS:i:87', 'XN:i:0', 'XM:i:2', 'XO:i:0', 'XG:i:0', 'NM:i:2', 'MD:Z:1T30A16', 'YT:Z:UU\n']

def bam_lines(f):
    p = subprocess.Popen('samtools view -',stdin=f,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    header_lines = []
    for line in p.stdout:
        yield line.decode('ascii')

def getBamReadPairs(bamfile1,bamfile2, skip_repeated_reads=False):
    bam1 = bam_lines(bamfile1)
    bam2 = bam_lines(bamfile2)
    line1= next(bam1).strip('\n').split() #split on white space. Results in 11 fields of mandatory SAM + variable number of additional tags.
    line2= next(bam2).strip('\n').split()
    while line1 and line2 and line1 !=[''] and line2 !=['']:
        #print(line1[0],line2[0],file=sys.stderr)
        assert line1[0] == line2[0]
        yield line1,line2
        previous_read1 = line1[0]
        previous_read2 = line2[0]
        if skip_repeated_reads:
            while line1 and line2 and line1 !=[''] and line2 !=[''] and line1[0] == previous_read1:
                line1= next(bam1).strip('\n').split()
            while line1 and line2 and line1 !=[''] and line2 !=[''] and line2[0] == previous_read2:
                line2= next(bam2).strip('\n').split()
        else:
            line1= next(bam1).strip('\n').split()
            line2= next(bam2).strip('\n').split()
    pass

def getReadPairs(sam1,sam2, skip_repeated_reads=False):
    line1= sam1.readline().strip('\n').split() #split on white space. Results in 11 fields of mandatory SAM + variable number of additional tags.
    line2= sam2.readline().strip('\n').split()
    while line1 and line2 and line1 !=[''] and line2 !=['']:
        #print(line1[0],line2[0],file=sys.stderr)
        assert line1[0] == line2[0]
        yield line1,line2
        previous_read1 = line1[0]
        previous_read2 = line2[0]
        if skip_repeated_reads:
            while line1 and line2 and line1 !=[''] and line2 !=[''] and line1[0] == previous_read1:
                line1= sam1.readline().strip('\n').split()
            while line1 and line2 and line1 !=[''] and line2 !=[''] and line2[0] == previous_read2:
                line2= sam2.readline().strip('\n').split()
        else:
            line1= sam1.readline().strip('\n').split()
            line2= sam2.readline().strip('\n').split()
            
    pass

def process_headers(file1,file2, primary_specific=sys.stdout, secondary_specific=None, primary_multi=None, secondary_multi=None, unassigned=None, unresolved=None, bam=False):
    if bam:
        samheader1 = "\n".join(get_bam_header(file1))
        samheader2 = "\n".join(get_bam_header(file2))
    else:
        samheader1 = "\n".join(get_sam_header(file1))
        samheader2 = "\n".join(get_sam_header(file2))
    print(samheader1, file=primary_specific)
    if secondary_specific:
        print(samheader2, file=secondary_specific)
    if primary_multi:
        print(samheader1, file=primary_multi)
    if secondary_multi:
        print(samheader2, file=secondary_multi)
    if unassigned:
        print(samheader1, file=unassigned)
    if unresolved: #This will not be the correct header - Look into merging header 1 and 2
        print(samheader1, file=unresolved)
    pass

def main_single_end(readpairs, primary_specific=sys.stdout, secondary_specific=None, primary_multi=None, secondary_multi=None, unassigned=None, unresolved=None):
    #assume that reads occur only once and are in the same order in both files
    #TO DO: add a test for these conditions
    for line1,line2 in readpairs:
        assert line1[0] == line2[0]
        AS1 = [int(x.split(':')[-1]) for x in line1[11:] if 'AS' in x]
        XS1 = [int(x.split(':')[-1]) for x in line1[11:] if 'XS' in x]
        AS2 = [int(x.split(':')[-1]) for x in line2[11:] if 'AS' in x]
        XS2 = [int(x.split(':')[-1]) for x in line2[11:] if 'XS' in x]
        if AS1 and (not AS2 or AS1 > AS2): #maps in primary better than secondary
            if not XS1 or AS1 > XS1:       #maps uniquely in primary better than secondary
                print('\t'.join(line1),file=primary_specific) #12 tab separated fields, all additional tags space sep.
            elif primary_multi:            #multimaps in primary better than secondary
                print('\t'.join(line1),file=primary_multi) 
        elif not AS2 and not AS1:          #does not map in either
            if unassigned: print('\t'.join(line1),file=unassigned)
        elif AS1 == AS2:                   #maps equally well in both
            if unresolved:
                print('\t'.join(line1),file=unresolved) #does this need to be two unresolved mapping files?
                print('\t'.join(line2),file=unresolved)
        elif AS2 and ((not AS1) or AS2 > AS1): #maps in secondary better than primary
            if (not XS2) or AS2 > XS2:
                if secondary_specific:
                    print('\t'.join(line2),file=secondary_specific) #maps uniquely in secondary better than primary
            elif secondary_multi:
                    print('\t'.join(line2),file=secondary_multi) #multimaps in secondary better than primary
    pass

def main_paired_end(readpairs, primary_specific=sys.stdout, secondary_specific=None, primary_multi=None, secondary_multi=None, unassigned=None, unresolved=None):
    #assume that paired end reads are sequential in the bam file, occur only once, and are in the same order in both files
    #TO DO: add a test for these conditions
    
    #Potential states:
    #   pair maps uniquely to primary = primary_specific
    #   pair maps uniquely to secondary = secondary_specific
    #   #Note: do not currently require mapping as a pair - just need both ends to single end map
    #   one read multimaps to primary and second read maps uniquely to primary = primary_specific
    #   one read multimaps to secondary and second read maps uniquely to secondary = secondary_specific
    #   both reads multimap to primary = primary_multi
    #   both reads multimap to secondary = secondary_multi
    #   one read does not map at all = unassigned (this is a conservative allocation - also puts virus/transgene boundary reads in this file)
    #   one read maps to primary, other read to secondary = unresolved

    previous_line1 = []
    previous_line2 = []
    for line1,line2 in readpairs:
        assert line1[0] == line2[0]
        
        #Test to see if line1 is pair of previous_line1 - if not skip ahead to next readpair
        if not previous_line1 or not line1 or not previous_line1[0] == line1[0]:
            previous_line1 = line1
            previous_line2 = line2
            continue
        
        #Get AS and XS tags from all four reads
        PAS1 = [int(x.split(':')[-1]) for x in previous_line1[11:] if 'AS' in x]
        PXS1 = [int(x.split(':')[-1]) for x in previous_line1[11:] if 'XS' in x]
        PAS2 = [int(x.split(':')[-1]) for x in previous_line2[11:] if 'AS' in x]
        PXS2 = [int(x.split(':')[-1]) for x in previous_line2[11:] if 'XS' in x]
        AS1 = [int(x.split(':')[-1]) for x in line1[11:] if 'AS' in x]
        XS1 = [int(x.split(':')[-1]) for x in line1[11:] if 'XS' in x]
        AS2 = [int(x.split(':')[-1]) for x in line2[11:] if 'AS' in x]
        XS2 = [int(x.split(':')[-1]) for x in line2[11:] if 'XS' in x]
        
        if AS1 and (not AS2 or AS1 > AS2): #maps in primary better than secondary
            if not XS1 or AS1 > XS1:       #maps uniquely in primary better than secondary
                print('\t'.join(previous_line1[:11])+'\t'+' '.join(previous_line1[11:]),file=primary_specific) 
                print('\t'.join(line1),file=primary_specific) 
            elif primary_multi:            #multimaps in primary better than secondary
                print('\t'.join(previous_line1[:11])+'\t'+' '.join(previous_line1[11:]),file=primary_multi) 
                print('\t'.join(line1),file=primary_multi) 
        elif not AS2 and not AS1:          #does not map in either
            if unassigned: print('\t'.join(line1),file=unassigned)
        elif AS1 == AS2:                   #maps equally well in both
            if unresolved:
                print('\t'.join(previous_line1[:11])+'\t'+' '.join(previous_line1[11:]),file=unresolved) #does this need to be two unresolved mapping files?
                print('\t'.join(line1),file=unresolved)
                print('\t'.join(previous_line2[:11])+'\t'+' '.join(previous_line2[11:]),file=unresolved)
                print('\t'.join(line2),file=unresolved)
        elif AS2 and ((not AS1) or AS2 > AS1): #maps in secondary better than primary
            if (not XS2) or AS2 > XS2:
                if secondary_specific:
                    print('\t'.join(previous_line2[:11])+'\t'+' '.join(previous_line2[11:]),file=secondary_specific) #maps uniquely in secondary better than primary
                    print('\t'.join(line2),file=secondary_specific) #maps uniquely in secondary better than primary
            elif secondary_multi:
                    print('\t'.join(previous_line2[:11])+'\t'+' '.join(previous_line2[11:]),file=secondary_multi) #multimaps in secondary better than primary
                    print('\t'.join(line2),file=secondary_multi) #multimaps in secondary better than primary
        previous_line1 = line1
        previous_line2 = line2
    pass

def command_line_interface(*args,**kw):
    parser = argparse.ArgumentParser(prog = "xenomapper",
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description=textwrap.dedent("""\
                    A script for parsing pairs of sam files and returning sam files
                    containing only reads where no better mapping exist in other files.
                    Used for filtering reads where multiple species may contribute 
                    (eg human tissue xenografted into mouse, pathogen growing on plant).
                    """),
                    epilog = textwrap.dedent("""\
                    To use bam files in a bash shell use process subtitution:
                        xenomapper --primary_specific >(samtools view -bS - > outfilename.bam) \\
                                   --primary_sam      <(samtools view primary.bam) \\
                                   --secondary_sam    <(samtools view secondary.bam)
                    
                    This program is distributed in the hope that it will be useful,
                    but WITHOUT ANY WARRANTY; without even the implied warranty of
                    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n
                    
                    """),
                    )
    parser.add_argument('--primary_sam',
                        type=argparse.FileType('rt'),
                        default=None,
                        help='a SAM format Bowtie2 mapping output file corrisponding to the primary species of interest')
    parser.add_argument('--secondary_sam',
                        type=argparse.FileType('rt'),
                        default=None,
                        help='a SAM format Bowtie2 mapping output file corrisponding to the secondary or contaminating species')
    parser.add_argument('--primary_bam',
                        type=argparse.FileType('rb'),
                        default=None,
                        help='a BAM format Bowtie2 mapping output file corrisponding to the primary species of interest')
    parser.add_argument('--secondary_bam',
                        type=argparse.FileType('rb'),
                        default=None,
                        help='a BAM format Bowtie2 mapping output file corrisponding to the secondary or contaminating species')
    parser.add_argument('--primary_specific',
                        type=argparse.FileType('wt'),
                        default=sys.stdout,
                        help='name for SAM format output file for reads mapping to a specific location in the primary species')
    parser.add_argument('--secondary_specific',
                        type=argparse.FileType('wt'),
                        default=None,
                        help='name for SAM format output file for reads mapping to a specific location in the secondary species')
    parser.add_argument('--primary_multi',
                        type=argparse.FileType('wt'),
                        default=None,
                        help='name for SAM format output file for reads multi mapping in the primary species')
    parser.add_argument('--secondary_multi',
                        type=argparse.FileType('wt'),
                        default=None,
                        help='name for SAM format output file for reads multi mapping in the secondary species')
    parser.add_argument('--unassigned',
                        type=argparse.FileType('wt'),
                        default=None,
                        help='name for SAM format output file for unassigned (non-mapping) reads')
    parser.add_argument('--unresolved',
                        type=argparse.FileType('wt'),
                        default=None,
                        help='name for SAM format output file for unresolved (maps equally well in both species) reads')
    parser.add_argument('--paired',
                        action='store_true',
                        help='the SAM files consist of paired reads with forward and reverse reads occuring once and interlaced')
    parser.add_argument('--version',
                        action='store_true',
                        help='print version information and exit')
    args = parser.parse_args()
    if args.version:
        print(__version__)
        sys.exit()
    if (not args.primary_sam or not args.secondary_sam) and \
        (not args.primary_bam or not args.secondary_bam):
        print('ERROR: You must provide --primary_sam and --secondary_sam\n or --primary_bam and --secondary_bam\n')
        parser.print_help()
        sys.exit(1)
    return args
    

def main():
    args = command_line_interface()
    #print(args, file=sys.stderr)
    
    if args.primary_sam:
        process_headers(args.primary_sam,args.secondary_sam,
                            primary_specific=args.primary_specific,
                            secondary_specific=args.secondary_specific,
                            primary_multi=args.primary_multi,
                            secondary_multi=args.secondary_multi,
                            unassigned=args.unassigned,
                            unresolved=args.unresolved)
                        
        readpairs = getReadPairs(args.primary_sam,args.secondary_sam)
    else:
        process_headers(args.primary_bam,args.secondary_bam,
                            primary_specific=args.primary_specific,
                            secondary_specific=args.secondary_specific,
                            primary_multi=args.primary_multi,
                            secondary_multi=args.secondary_multi,
                            unassigned=args.unassigned,
                            unresolved=args.unresolved,
                            bam=True)
                        
        readpairs = getBamReadPairs(args.primary_bam,args.secondary_bam)
        
    
    if args.paired:
        main_paired_end(readpairs,
                        primary_specific=args.primary_specific,
                        secondary_specific=args.secondary_specific,
                        primary_multi=args.primary_multi,
                        secondary_multi=args.secondary_multi,
                        unassigned=args.unassigned,
                        unresolved=args.unresolved)
        
    else:
        main_single_end(readpairs,
                        primary_specific=args.primary_specific,
                        secondary_specific=args.secondary_specific,
                        primary_multi=args.primary_multi,
                        secondary_multi=args.secondary_multi,
                        unassigned=args.unassigned,
                        unresolved=args.unresolved)
    pass


if __name__ == '__main__':
    main()
    
