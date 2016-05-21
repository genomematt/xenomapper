#!/usr/bin/env python3
# encoding: utf-8
"""
xenomapper.py

A script for parsing pairs of sam files and returning sam files containing only reads where no better mapping exist in other files.
Used for filtering reads where multiple species may contribute (eg human tissue xenografted into mouse).

Created by Matthew Wakefield.
Copyright (c) 2011-2016  Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne. All rights reserved.

   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


"""
import sys
import os
import argparse, textwrap
import subprocess
import re
from collections import Counter
from copy import copy

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2016 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production/Stable"

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

def get_bam_header(bamfile): #pragma: no cover #not tested due to need for samtools
    p = subprocess.Popen('samtools view -H -',stdin=bamfile,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    header = []
    for line in p.stdout:
        header.append(line.decode('ascii'))
    bamfile.seek(0) #reset to start of file for next samtools call
    return [x.strip('\n') for x in header]

def bam_lines(f): #pragma: no cover #not tested due to need for samtools
    """Use samtools in a subprocess to yield lines of sam from a bam file
        Arguments: a file or file like object in bam format
        Yields:    ascii sam file lines
    """
    p = subprocess.Popen('samtools view -',stdin=f,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    header_lines = []
    for line in p.stdout:
        yield line.decode('ascii')

def getBamReadPairs(bamfile1,bamfile2, skip_repeated_reads=False): #pragma: no cover #not tested due to need for samtools
    """Process two bamfiles to yield the equivalent line from each file
        Arguments: 
        bamfile1, bamfile2  - file or file like objects in binary bam format
                              containing the same reads in the same order
                              mapped in two different species
        Yields:    a tuple of lists of sam fields split on white space
    """
    bam1 = bam_lines(bamfile1)
    bam2 = bam_lines(bamfile2)
    line1= next(bam1).strip('\n').split() #split on white space. Results in 11 fields of mandatory SAM + variable number of additional tags.
    line2= next(bam2).strip('\n').split()
    while line1 and line2 and line1 !=[''] and line2 !=['']:
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
    """Process two sam files to yield the equivalent line from each file
        Arguments: 
        sam1, sam2  - file or file like objects in ascii sam format
                      containing the same reads in the same order
                      mapped in two different species
        Yields:    a tuple of lists of sam fields split on white space
    """
    line1= sam1.readline().strip('\n').split() #split on white space. Results in 11 fields of mandatory SAM + variable number of additional tags.
    line2= sam2.readline().strip('\n').split()
    while line1 and line2 and line1 !=[''] and line2 !=['']:
        assert line1[0] == line2[0]
        yield line1,line2
        previous_read1 = line1[0]
        previous_read2 = line2[0]
        if skip_repeated_reads: #pragma: no cover
            while line1 and line2 and line1 !=[''] and line2 !=[''] and line1[0] == previous_read1:
                line1= sam1.readline().strip('\n').split()
            while line1 and line2 and line1 !=[''] and line2 !=[''] and line2[0] == previous_read2:
                line2= sam2.readline().strip('\n').split()
        else:
            line1= sam1.readline().strip('\n').split()
            line2= sam2.readline().strip('\n').split()
    pass

def add_pg_tag(sam_header_list,comment=None):
    new_header = copy(sam_header_list)
    if not [x[0] for x in new_header] == ['@',]*len(new_header):
        raise ValueError('Incorrect SAM header format :\n{0}'.format('\n'.join(new_header))) #pragma: no cover
    if new_header[-1][:3] == '@PG':
        PP = 'PP'+[x for x in new_header[-1].split() if x[:2] == 'ID'][0][2:]+'\t'
    else:
        PP = ''
    new_header.append('@PG\tID:Xenomapper\tPN:Xenomapper\t'+PP+'VN:{0}'.format(__version__))
    if comment:
        new_header.append('@CO\t'+comment)
    return new_header

def process_headers(file1,file2, primary_specific=sys.stdout, secondary_specific=None, primary_multi=None, secondary_multi=None, unassigned=None, unresolved=None, bam=False):
    """Process headers from two sam or bam files and write appropriate
    header information to the correct output files
        Arguments: 
        file1, file2  - file or file like objects in binary bam format
                        or ascii sam format
        primary_specific, secondary_specific, primary_multi,
        secondary_multi, unassigned, unresolved
                      - ascii file or file like objects for outputs
        bam           - Boolean flag indicating file1 & file2 are
                        in binary bam format.  Default = False
    """
    if bam: #pragma: no cover
        samheader1 = get_bam_header(file1)
        samheader2 = get_bam_header(file2)
    else:
        samheader1 = get_sam_header(file1)
        samheader2 = get_sam_header(file2)
    print('\n'.join(add_pg_tag(samheader1,
                    comment='species specific reads'
                    )), file=primary_specific)
    if secondary_specific:
        print('\n'.join(add_pg_tag(samheader2,
                    comment='species specific reads'
                    )), file=secondary_specific)
    if primary_multi:
        print('\n'.join(add_pg_tag(samheader1,
                    comment='species specific multimapping reads'
                    )), file=primary_multi)
    if secondary_multi:
        print('\n'.join(add_pg_tag(samheader2,
                    comment='species specific multimapping reads'
                    )), file=secondary_multi)
    if unassigned:
        print('\n'.join(add_pg_tag(samheader1,
                    comment='reads that could not be assigned'
                    )), file=unassigned)
    if unresolved: #This will not be the correct header - Look into merging header 1 and 2
        print('\n'.join(add_pg_tag(samheader1,
                    comment='reads that could not be resolved'
                    )), file=unresolved)
    pass

def get_tag(sam_line,tag='AS'):
    """Return the value of a SAM tag field
    Arguments:
        sam_line  - list of elements from a SAM file line
        tag       - the name of the optional tag to be returned
                    Only suitable for numeric tags eg AS or XS
    Returns
        tag_value - the value of the SAM tag converted to a float
                    or -inf if tag is not present.
    """
    tag_list = [x for x in sam_line[11:] if tag in x]
    if not tag_list:
        return float('-inf') #this will always be worse than any bowtie score
    if len(tag_list) > 1:
        raise ValueError('SAM line has multiple values of {0}: {1}'.format(tag,sam_line)) #pragma: no cover
    return float(tag_list[0].split(':')[-1])

def get_tag_with_ZS_as_XS(sam_line,tag='AS'):
    """Return the value of a SAM tag field
    Arguments:
        sam_line  - list of elements from a SAM file line
        tag       - the name of the optional tag to be returned
                    Only suitable for numeric tags eg AS or XS
                    Requests for XS will return the value for ZS
    Returns
        tag_value - the value of the SAM tag converted to a float
                    or -inf if tag is not present.
    """
    if tag == 'XS':
        tag = 'ZS' 
    return get_tag(sam_line,tag)

#def alternative_get_tag(sam_line,tag='AS'):
#    """Dummy code example for overriding get_tag to support
#    aligners that produce a useful score but store it in an
#    alternive tag. Our hypothetical aligner uses ZS for
#    match score and ZM for next highest hit.
#    """
#    if tag == 'AS':
#        other_tag = 'ZS'
#    elif tag == 'XS':
#        other_tag = 'ZM'
#    else:
#        return get_tag(sam_line,tag)
#    tag_list = [x for x in sam_line[11:] if other_tag in x]
#    if not tag_list:
#        return float('-inf') #this will always be worse than any bowtie score
#    if len(tag_list) > 1:
#        raise ValueError('SAM line has multiple values of {0}: {1}'.format(other_tag,sam_line))
#    return float(tag_list[0].split(':')[-1])


def get_cigarbased_AS_tag(sam_line,tag='AS'):
    """Return calculated values for a score equivalent to a
    rescaled AS tag using the cigar line and NM tags from the
    SAM entry.
    Score is rescaled such that a perfect match for the full
    length of the read acchieves a score of 0.0.
    Penalties are: mismatch -6, gap open -5, gap extend -3,
                   softclipped -2 (equiv to +2 match rescaled)
    Arguments:
        sam_line  - list of elements from a SAM file line
        tag       - the name of the optional tag to be emulated
                    Only AS tags will be calculated
    Returns
        tag_value - the calculated value of the AS score as a
                    float or the return value of get_tag for
                    all other values of tag
    """
    if tag != 'AS':
        return get_tag(sam_line,tag)
    NM = [x for x in sam_line[11:] if 'NM' in x]
    if not NM:
        return float('-inf') #either a multimapper or unmapped
    mismatches = int(NM[0].split(':')[-1])
    cigar = re.findall(r'([0-9]+)([MIDNSHPX=])',sam_line[5]) 
    deletions = [int(x[0]) for x in cigar if x[1] == 'D']
    insertions = [int(x[0]) for x in cigar if x[1] == 'I']
    softclips = [int(x[0]) for x in cigar if x[1] == 'S']
    score = (-6 * mismatches) + (-5 * (len(insertions) + len(deletions))) + (-3 * (sum(insertions) + sum(deletions))) + (-2 * sum(softclips))
    return score

def get_mapping_state(AS1,XS1,AS2,XS2, min_score=float('-inf')):
    """Determine the mapping state based on scores in each species.
    Scores can be negative but better matches must have higher scores
    Arguments:
        AS1  - float or interger score of best match in primary species
        XS1  - float or interger score of other match in primary species
        AS2  - float or interger score of best match in secondary species
        XS2  - float or interger score of other match in secondary species
        min_score - the score that matches must exceed in order to be
                    considered valid matches. Note scores equalling this
                    value will also be considered not to match.
                    Default = -inf
    Returns
        state - a string of 'primary_specific', 'secondary_specific',
                'primary_multi', 'secondary_multi',
                'unresolved', or 'unassigned' indicating match state.
    """
    if AS1 <= min_score and  AS2 <= min_score:  #low quality mapping in both
        return 'unassigned'
    elif AS1 > min_score and (AS2 <= min_score or AS1 > AS2): #maps in primary better than secondary
        if not XS1 or AS1 > XS1:       #maps uniquely in primary better than secondary
            return 'primary_specific'
        else:            #multimaps in primary better than secondary
            return 'primary_multi' 
    elif AS1 == AS2:                   #maps equally well in both
        return 'unresolved'
    elif AS2 > min_score and (AS1 <= min_score or AS2 > AS1): #maps in secondary better than primary
        if (not XS2) or AS2 > XS2:
            return 'secondary_specific'
        else:
            return 'secondary_multi' #multimaps in secondary better than primary
    else: raise RuntimeError('Error in processing logic with values {0} '.format((AS1,XS1,AS2,XS2))) # pragma: no cover

def main_single_end(readpairs,
                    primary_specific=sys.stdout,
                    secondary_specific=None,
                    primary_multi=None,
                    secondary_multi=None,
                    unassigned=None,
                    unresolved=None,
                    min_score=float('-inf'),
                    tag_func=get_tag):
    """Main loop for processing single end read files
    Arguments:
        readpairs - an iterable of tuples of lists of sam fields
        primary_specific, secondary_specific, primary_multi,
        secondary_multi, unassigned, unresolved
                  - ascii file or file like objects for outputs
        min_score - the score that matches must exceed in order to be
                    considered valid matches. Note scores equalling this
                    value will also be considered not to match.
                    Default = -inf
        tag_func  - a function that takes a list of sam fields and a
                    tag identifier (at least 'AS' and 'XS')
                    returns a numeric value for that tag
    Returns:
        category_counts - a dictionary keyed by category containing
                    occurance counts
    """
    #assume that reads occur only once and are in the same order in both files
    
    category_counts = Counter()
    
    for line1,line2 in readpairs:
        assert line1[0] == line2[0]
        AS1 = tag_func(line1, tag='AS')
        XS1 = tag_func(line1, tag='XS')
        AS2 = tag_func(line2, tag='AS')
        XS2 = tag_func(line2, tag='XS')
        
        state = get_mapping_state(AS1,XS1,AS2,XS2,min_score)
        
        category_counts[state] += 1
        
        if state == 'primary_specific':
            if primary_specific:
                print('\t'.join(line1),file=primary_specific)
        elif state == 'secondary_specific':
            if secondary_specific:
                print('\t'.join(line2),file=secondary_specific)
        elif state == 'primary_multi':
            if primary_multi:
                print('\t'.join(line1),file=primary_multi)
        elif state == 'secondary_multi':
            if secondary_multi:
                print('\t'.join(line2),file=secondary_multi)
        elif state == 'unassigned':
            if unassigned:
                print('\t'.join(line1),file=unassigned)
        elif state == 'unresolved':
            if unresolved:
                print('\t'.join(line1),file=unresolved)
                print('\t'.join(line2),file=unresolved)
        else: raise RuntimeError('Unexpected state {0} '.format(state)) # pragma: no cover
    return category_counts

def main_paired_end(readpairs,
                    primary_specific=sys.stdout,
                    secondary_specific=None,
                    primary_multi=None,
                    secondary_multi=None,
                    unassigned=None,
                    unresolved=None,
                    min_score=float('-inf'),
                    tag_func=get_tag):
    """Liberal main loop for processing paired end read files.
    Discordant reads will be assigned to the highest
    priority category in the order: primary_specific,
    secondary_specific, primary_multi, secondary_multi,
    unresolved, unassigned
    
    This will rescue unresolved or unassigned reads where the second read 
    can be assigned, and deems primary/secondary discordant reads as
    primary.
    
    Paired end reads must be sequential in the sam or bam file,
    occur only once and be in the same order in both files.
    
    Arguments:
        readpairs - an iterable of tuples of lists of sam fields
        primary_specific, secondary_specific, primary_multi,
        secondary_multi, unassigned, unresolved
                  - ascii file or file like objects for outputs
        min_score - the score that matches must exceed in order to be
                    considered valid matches. Note scores equalling this
                    value will also be considered not to match.
                    Default = -inf
        tag_func  - a function that takes a list of sam fields and a
                    tag identifier (at least 'AS' and 'XS')
                    returns a numeric value for that tag
    Returns:
        category_counts - a dictionary keyed by a tuple of forward
                    and reverse read category containing
                    occurance counts
    """
    
    category_counts = Counter()
    
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
        PAS1 = tag_func(previous_line1, tag='AS')
        PXS1 = tag_func(previous_line1, tag='XS')
        PAS2 = tag_func(previous_line2, tag='AS')
        PXS2 = tag_func(previous_line2, tag='XS')
        AS1 = tag_func(line1, tag='AS')
        XS1 = tag_func(line1, tag='XS')
        AS2 = tag_func(line2, tag='AS')
        XS2 = tag_func(line2, tag='XS')
        
        forward_state = get_mapping_state(PAS1,PXS1,PAS2,PXS2,min_score)
        reverse_state = get_mapping_state(AS1,XS1,AS2,XS2,min_score)
        
        category_counts[(forward_state,reverse_state)] += 1
        
        
        if forward_state == 'primary_specific' or reverse_state == 'primary_specific':
            if primary_specific:
                print('\t'.join(previous_line1),file=primary_specific) 
                print('\t'.join(line1),file=primary_specific)
        elif forward_state == 'secondary_specific' or reverse_state == 'secondary_specific':
            if secondary_specific:
                print('\t'.join(previous_line2),file=secondary_specific) 
                print('\t'.join(line2),file=secondary_specific)
        elif forward_state == 'primary_multi' or reverse_state == 'primary_multi':
            if primary_multi:
                print('\t'.join(previous_line1),file=primary_multi) 
                print('\t'.join(line1),file=primary_multi)
        elif forward_state == 'secondary_multi' or reverse_state == 'secondary_multi':
            if secondary_multi:
                print('\t'.join(previous_line2),file=secondary_multi) 
                print('\t'.join(line2),file=secondary_multi)
        elif forward_state == 'unresloved' or reverse_state == 'unresolved':
            if unresolved:
                print('\t'.join(previous_line1),file=unresolved) 
                print('\t'.join(line1),file=unresolved)
                print('\t'.join(previous_line2),file=unresolved) 
                print('\t'.join(line2),file=unresolved)
        elif forward_state == 'unassigned' or reverse_state == 'unassigned':
            if unassigned:
                print('\t'.join(previous_line1),file=unassigned) 
                print('\t'.join(line1),file=unassigned)
        else: raise RuntimeError('Unexpected states forward:{0} reverse:{1}'.format(forward_state,reverse_state)) # pragma: no cover
        
        previous_line1 = line1
        previous_line2 = line2
        
    return category_counts

def conservative_main_paired_end(readpairs,
                    primary_specific=sys.stdout,
                    secondary_specific=None,
                    primary_multi=None,
                    secondary_multi=None,
                    unassigned=None,
                    unresolved=None,
                    min_score=float('-inf'),
                    tag_func=get_tag):
    """Main loop for conservative processing of paired end read files.
    Read pairs where either read is unassigned will be deemed unassigned.
    This places features such as transgene boundaries in the unassigned file.
    Reads with discordant species will be unresolved. 
    Concordant reads will be assigned to the highest priority category in the 
    order: primary_specific, secondary_specific, primary_multi, 
    secondary_multi, unresolved, unassigned
    
    Paired end reads must be sequential in the sam or bam file, occur only once
    and be in the same order in both files.
    
    Arguments:
        readpairs - an iterable of tuples of lists of sam fields
        primary_specific, secondary_specific, primary_multi,
        secondary_multi, unassigned, unresolved
                  - ascii file or file like objects for outputs
        min_score - the score that matches must exceed in order to be
                    considered valid matches. Note scores equalling this
                    value will also be considered not to match.
                    Default = -inf
        tag_func  - a function that takes a list of sam fields and a
                    tag identifier (at least 'AS' and 'XS')
                    returns a numeric value for that tag
    Returns:
        category_counts - a dictionary keyed by a tuple of forward
                    and reverse read category containing
                    occurance counts
    """
    
    category_counts = Counter()
    
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
        PAS1 = tag_func(previous_line1, tag='AS')
        PXS1 = tag_func(previous_line1, tag='XS')
        PAS2 = tag_func(previous_line2, tag='AS')
        PXS2 = tag_func(previous_line2, tag='XS')
        AS1 = tag_func(line1, tag='AS')
        XS1 = tag_func(line1, tag='XS')
        AS2 = tag_func(line2, tag='AS')
        XS2 = tag_func(line2, tag='XS')
        
        forward_state = get_mapping_state(PAS1,PXS1,PAS2,PXS2,min_score)
        reverse_state = get_mapping_state(AS1,XS1,AS2,XS2,min_score)
        
        category_counts[(forward_state,reverse_state)] += 1
        if forward_state == 'unassigned' or reverse_state == 'unassigned':
            if unassigned:
                print('\t'.join(previous_line1),file=unassigned) 
                print('\t'.join(line1),file=unassigned)
        elif forward_state == 'unresloved' or reverse_state == 'unresolved' \
            or (forward_state in ['primary_specific','primary_multi'] and \
                reverse_state in ['secondary_specific','secondary_multi']) \
            or (forward_state in ['secondary_specific','secondary_multi'] and \
                reverse_state in ['primary_specific','primary_multi']):
            if unresolved:
                print('\t'.join(previous_line1),file=unresolved) 
                print('\t'.join(line1),file=unresolved)
                print('\t'.join(previous_line2),file=unresolved) 
                print('\t'.join(line2),file=unresolved)
        elif forward_state == 'primary_specific' or reverse_state == 'primary_specific':
            if primary_specific:
                print('\t'.join(previous_line1),file=primary_specific) 
                print('\t'.join(line1),file=primary_specific)
        elif forward_state == 'secondary_specific' or reverse_state == 'secondary_specific':
            if secondary_specific:
                print('\t'.join(previous_line2),file=secondary_specific) 
                print('\t'.join(line2),file=secondary_specific)
        elif forward_state == 'primary_multi' or reverse_state == 'primary_multi':
            if primary_multi:
                print('\t'.join(previous_line1),file=primary_multi) 
                print('\t'.join(line1),file=primary_multi)
        elif forward_state == 'secondary_multi' or reverse_state == 'secondary_multi':
            if secondary_multi:
                print('\t'.join(previous_line2),file=secondary_multi) 
                print('\t'.join(line2),file=secondary_multi)
        else: raise RuntimeError('Unexpected states forward:{0} reverse:{1}'.format(forward_state,reverse_state)) # pragma: no cover
        
        previous_line1 = line1
        previous_line2 = line2
        
    return category_counts

def output_summary(category_counts, outfile=sys.stderr):
    print('-'*80, file=outfile)
    print('Read Count Category Summary\n', file=outfile)
    print('|       {0:45s}|     {1:10s}  |'.format('Category','Count'), file=outfile)
    print('|:','-'*50,':|:','-'*15,':|',sep='', file=outfile)
    for category in sorted(category_counts):
        print('|  {0:50s}|{1:15d}  |'.format(str(category),category_counts[category]), file=outfile)
    print(file=outfile)
    pass

def command_line_interface(*args,**kw): #pragma: no cover
    parser = argparse.ArgumentParser(prog = "xenomapper",
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description=textwrap.dedent("""\
                    A script for parsing pairs of sam files and returning sam files
                    containing only reads where no better mapping exist in other files.
                    Used for filtering reads where multiple species may contribute 
                    (eg human tissue xenografted into mouse, pathogen growing on plant).
                    
                    Files should contain an AS and XS score and better matches must have
                    a higher alignment score (but can be negative).
                    Reads must be in the same order in both species.
                    
                    In practice this is best acchieved by using Bowtie2 in --local mode.
                    If the -p option is used you must also use --reorder.
                    
                    Limited support is provided for aligners that do not produce AS and XS
                    score tags via the --cigar_score option.
                    
                    All input files must be seekable
                    (ie not a FIFO, process substitution or pipe)'
                    """),
                    epilog = textwrap.dedent("""\
                    To output bam files in a bash shell use process subtitution:
                        xenomapper --primary_specific >(samtools view -bS - > outfilename.bam) 
                    
                    This program is distributed in the hope that it will be useful,
                    but WITHOUT ANY WARRANTY; without even the implied warranty of
                    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n
                    
                    """),
                    )
    parser.add_argument('--primary_sam',
                        type=argparse.FileType('rt'),
                        default=None,
                        help='a SAM format Bowtie2 mapping output file corresponding to the primary species of interest')
    parser.add_argument('--secondary_sam',
                        type=argparse.FileType('rt'),
                        default=None,
                        help='a SAM format Bowtie2 mapping output file corresponding to the secondary or contaminating species')
    parser.add_argument('--primary_bam',
                        type=argparse.FileType('rb'),
                        default=None,
                        help='a BAM format Bowtie2 mapping output file corresponding to the primary species of interest')
    parser.add_argument('--secondary_bam',
                        type=argparse.FileType('rb'),
                        default=None,
                        help='a BAM format Bowtie2 mapping output file corresponding to the secondary or contaminating species')
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
    parser.add_argument('--conservative',
                        action='store_true',
                        help='conservatively allocate paired end reads with discordant category allocations. \
                              Only pairs that are both specific, or specific and multi will be allocated as specific. \
                              Pairs that are discordant for species will be deemed unresolved.  Pairs where any read \
                              is unassigned will be deemed unassigned.')
    parser.add_argument('--min_score',
                        type=float,
                        default=float('-inf'),
                        help='the minimum mapping score.  Reads with scores less than or equal to min_score will be considered unassigned. \
                              Values should be chosen based on the mapping program and read length')
    parser.add_argument('--cigar_scores',
                        action='store_true',
                        help='Use the cigar line and the NM tag to calculate a score. For aligners that do not support the AS tag. \
                              No determination of multimapping state will be done.  Reads that are unique in one species and multimap \
                              in the other species may be misassigned as no score can be calculated in the multimapping species. \
                              Score is -6 * mismatches + -5 * indel open + -3 * indel extend + -2 * softclip. \
                              Treatment of multimappers will vary with aligner.  If multimappers are assigned a cigar line they \
                              will be treated as species specific, otherwise as unassigned.')
    parser.add_argument('--use_zs',
                        action='store_true',
                        help='Use the value of the ZS tag in place of XS for determining the mapping score of the next best \
                              alignment.  Used with HISAT as the XS:A tag is conventionally used for strand in spliced mappers.')
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
    

def main(): #pragma: no cover
    args = command_line_interface()
    
    if args.cigar_scores:
        tag_func = get_cigarbased_AS_tag
    elif args.use_zs:
        tag_func = get_tag_with_ZS_as_XS
    else:
        tag_func = get_tag
    
    skip_repeated = False if args.paired else True
    
    if args.primary_sam:
        process_headers(args.primary_sam,args.secondary_sam,
                            primary_specific=args.primary_specific,
                            secondary_specific=args.secondary_specific,
                            primary_multi=args.primary_multi,
                            secondary_multi=args.secondary_multi,
                            unassigned=args.unassigned,
                            unresolved=args.unresolved)
                        
        readpairs = getReadPairs(args.primary_sam, args.secondary_sam, skip_repeated_reads=skip_repeated)
    else:
        process_headers(args.primary_bam,args.secondary_bam,
                            primary_specific=args.primary_specific,
                            secondary_specific=args.secondary_specific,
                            primary_multi=args.primary_multi,
                            secondary_multi=args.secondary_multi,
                            unassigned=args.unassigned,
                            unresolved=args.unresolved,
                            bam=True)
                        
        readpairs = getBamReadPairs(args.primary_bam, args.secondary_bam, skip_repeated_reads=skip_repeated)
        
    
    if args.paired:
        if args.conservative:
            paired_function = conservative_main_paired_end
        else:
            paired_function = main_paired_end
        category_counts = paired_function(readpairs,
                        primary_specific=args.primary_specific,
                        secondary_specific=args.secondary_specific,
                        primary_multi=args.primary_multi,
                        secondary_multi=args.secondary_multi,
                        unassigned=args.unassigned,
                        unresolved=args.unresolved,
                        min_score=args.min_score,
                        tag_func=tag_func)
        
    else:
        category_counts = main_single_end(readpairs,
                        primary_specific=args.primary_specific,
                        secondary_specific=args.secondary_specific,
                        primary_multi=args.primary_multi,
                        secondary_multi=args.secondary_multi,
                        unassigned=args.unassigned,
                        unresolved=args.unresolved,
                        min_score=args.min_score,
                        tag_func=tag_func)
    
    output_summary(category_counts=category_counts)
    pass


if __name__ == '__main__': #pragma: no cover
    main() 
    
