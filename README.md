[![Build Status](https://travis-ci.org/genomematt/xenomapper.svg?branch=master)](https://travis-ci.org/genomematt/xenomapper)
[![Coverage Status](https://coveralls.io/repos/genomematt/xenomapper/badge.svg)](https://coveralls.io/r/genomematt/xenomapper)

Xenomapper
==========

Xenomapper is a utility for post processing Bowtie2 (Langmead, 2012) mapped reads that have been aligned to a primary genome and a secondary genome and binning reads into species specific, multimapping in each species, unmapped and unassigned bins.  It can be used on single end or paired end sequencing data.  In paired end data it is assumed both reads come from the same species and evidence of sequence specificity for either read will be used to assign both reads.

Use cases include xenografts of human cancers in mouse and host pathogen interactions.

![Schematic of Xenomapper Use](/schematic.jpg "Schematic of Xenomapper Use")

Installation
============
Xenomapper requires python 3.4 or higher and is tested on linux and MacOS.  For bam file decoding samtools must be installed.


    git clone https://github.com/genomematt/xenomapper
    pip3 install xenomapper
    python3 -m xenomapper.tests.test_all

All users should upgrade to v0.5.0 or higher

Using Xenomapper
================

usage:

    xenomapper [-h]   [--primary_sam PRIMARY_SAM]
                      [--secondary_sam SECONDARY_SAM]
					  [--primary_bam PRIMARY_BAM]
                      [--secondary_bam SECONDARY_BAM]
                      [--primary_specific PRIMARY_SPECIFIC]
                      [--secondary_specific SECONDARY_SPECIFIC]
                      [--primary_multi PRIMARY_MULTI]
                      [--secondary_multi SECONDARY_MULTI]
                      [--unassigned UNASSIGNED]
					  [--unresolved UNRESOLVED]
                      [--paired]
					  [--min_score MIN_SCORE]
					  [--cigar_scores]
                      [--version]

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

	optional arguments:
	  -h, --help            show this help message and exit
	  --primary_sam PRIMARY_SAM
	                        a SAM format Bowtie2 mapping output file corrisponding
	                        to the primary species of interest
	  --secondary_sam SECONDARY_SAM
	                        a SAM format Bowtie2 mapping output file corrisponding
	                        to the secondary or contaminating species
	  --primary_bam PRIMARY_BAM
	                        a BAM format Bowtie2 mapping output file corrisponding
	                        to the primary species of interest
	  --secondary_bam SECONDARY_BAM
	                        a BAM format Bowtie2 mapping output file corrisponding
	                        to the secondary or contaminating species
	  --primary_specific PRIMARY_SPECIFIC
	                        name for SAM format output file for reads mapping to a
	                        specific location in the primary species
	  --secondary_specific SECONDARY_SPECIFIC
	                        name for SAM format output file for reads mapping to a
	                        specific location in the secondary species
	  --primary_multi PRIMARY_MULTI
	                        name for SAM format output file for reads multi
	                        mapping in the primary species
	  --secondary_multi SECONDARY_MULTI
	                        name for SAM format output file for reads multi
	                        mapping in the secondary species
	  --unassigned UNASSIGNED
	                        name for SAM format output file for unassigned (non-
	                        mapping) reads
	  --unresolved UNRESOLVED
	                        name for SAM format output file for unresolved (maps
	                        equally well in both species) reads
	  --paired              the SAM files consist of paired reads with forward and
	                        reverse reads occuring once and interlaced
	  --min_score MIN_SCORE
							the minimum mapping score.  Reads with scores less than
							or equal to min_score will be considered unassigned.
							Values should be chosen based on the mapping program 
							and read length
	  --cigar_scores        Use the cigar line and the NM tag to calculate a
	                        score. For aligners that do not support the AS tag. No
	                        determination of multimapping state will be done.
	                        Reads that are unique in one species and multimap in
	                        the other species may be misassigned as no score can
	                        be calculated in the multimapping species. Score is -6
	                        * mismatches + -5 * indel open + -3 * indel extend +
	                        -2 * softclip.
	  --version             print version information and exit


To output bam files in a bash shell use process substitution:


    xenomapper --primary_specific >(samtools view -bS - > outfilename.bam)


xenomappability
===============
xenomappability is a tool for creating mappability wiggle files that reflect the paired end and multi species nature of the final number more accurately than the commonly used single end mappability tracks.

This feature is computationally intensive for useful genomes.  In most cases you will want to segment into chromosomal or smaller regions and calculate on a cluster.


    xenomappability --fasta tests/data/test_from_EcoliK12DH10B.fasta --readlength 10 > tests/data/test_from_EcoliK12DH10B_10reads.fasta

    bowtie2-build tests/data/test_from_EcoliK12DH10B.fasta tests/data/test_from_EcoliK12DH10B
    bowtie2 -x tests/data/test_from_EcoliK12DH10B -f -U tests/data/test_from_EcoliK12DH10B_10reads.fasta -S tests/data/test_from_EcoliK12DH10B_10reads.sam

    xenomappability --mapped_test_data tests/data/test_from_EcoliK12DH10B_10reads.sam > tests/data/test_from_EcoliK12DH10B_10reads.wig
    xenomappability --single_end_wiggle tests/data/test_from_EcoliK12DH10B_10reads.wig --sam_for_sizes tests/data/paired_end_testdata_human.sam`
	
Citing Xenomapper
=================
Currently Xenomapper is unpublished, but this repository does have a DOI identifier for each release you can use to cite the code.  The DOI for the most current version is [![DOI](https://zenodo.org/badge/11450/genomematt/xenomapper.svg)](http://dx.doi.org/10.5281/zenodo.16677)

References
=================
Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359. http://bowtie-bio.sourceforge.net/bowtie2/

