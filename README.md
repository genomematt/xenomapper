[![Build Status](https://travis-ci.org/genomematt/xenomapper.svg?branch=master)](https://travis-ci.org/genomematt/xenomapper)
[![Coverage Status](https://coveralls.io/repos/genomematt/xenomapper/badge.svg)](https://coveralls.io/r/genomematt/xenomapper)
[![JOSS](http://joss.theoj.org/papers/7065e20e97ff4e44695b5e9fef7cdcd8/status.svg)](http://dx.doi.org/10.21105.joss.00018)
[![DOI](https://zenodo.org/badge/11450/genomematt/xenomapper.svg)](https://zenodo.org/badge/latestdoi/11450/genomematt/xenomapper)

Xenomapper
==========

Xenomapper is a utility for post processing mapped reads that have been aligned to a primary genome and a secondary genome and binning reads into species specific, multimapping in each species, unmapped and unassigned bins.  It can be used on single end or paired end sequencing data.  In paired end data evidence of sequence specificity for either read will be used to assign both reads.

Use cases include xenografts of human cancers and host pathogen interactions.

Xenomapper is most effective with mapped reads that include an XS or ZS score that gives the mapping score of the next best read.  These include Bowtie2 (Langmead, 2012) and HISAT (Kim, 2015). 

![Schematic of Xenomapper Use](/schematic.jpg "Schematic of Xenomapper Use")

Installation
============
Xenomapper requires python 3.3 or higher and is tested on linux and MacOS with CPython and pypy3.  For bam file decoding samtools must be installed.

Installing from the Python Package Index with pip is the easiest option:

    pip3 install xenomapper
    
Alternatively if you would like to install from the github repository

    git clone https://github.com/genomematt/xenomapper
    pip3 install --upgrade xenomapper
	
Although the repository tests by continuous integration with TravisCI its good practice to run the tests locally and check your install works correctly.  The tests are run with the following command:

    python3 -m xenomapper.tests.test_all

All users should upgrade to v0.5.0 or higher as earlier versions have known bugs.

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
	                        a SAM format Bowtie2 mapping output file corresponding
	                        to the primary species of interest
	  --secondary_sam SECONDARY_SAM
	                        a SAM format Bowtie2 mapping output file corresponding
	                        to the secondary or contaminating species
	  --primary_bam PRIMARY_BAM
	                        a BAM format Bowtie2 mapping output file corresponding
	                        to the primary species of interest
	  --secondary_bam SECONDARY_BAM
	                        a BAM format Bowtie2 mapping output file corresponding
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
	  --conservative        conservatively allocate paired end reads with
	                        discordant category allocations. Only pairs that are
	                        both specific, or specific and multi will be allocated
	                        as specific. Pairs that are discordant for species
	                        will be deemed unresolved. Pairs where any read is
	                        unassigned will be deemed unassigned.
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
	  --use_zs              Use the value of the ZS tag in place of XS for
	                        determining the mapping score of the next best
	                        alignment. Used with HISAT as the XS:A tag is
	                        conventionally used for strand in spliced mappers.
	  --version             print version information and exit


To output bam files in a bash shell use process substitution:


    xenomapper --primary_specific >(samtools view -bS - > outfilename.bam)


A worked example of using xenomapper can be found in [example_usage.ipynb](example_usage.ipynb)

xenomappability
===============
xenomappability is a tool for creating mappability wiggle files that reflect the paired end and multi species nature of the final number more accurately than the commonly used single end mappability tracks.

This feature is computationally intensive for useful genomes.  In most cases you will want to segment into chromosomal or smaller regions and calculate on a cluster.


    xenomappability --fasta tests/data/test_from_EcoliK12DH10B.fasta --readlength 10 > tests/data/test_from_EcoliK12DH10B_10reads.fasta

    bowtie2-build tests/data/test_from_EcoliK12DH10B.fasta tests/data/test_from_EcoliK12DH10B
    bowtie2 -x tests/data/test_from_EcoliK12DH10B -f -U tests/data/test_from_EcoliK12DH10B_10reads.fasta -S tests/data/test_from_EcoliK12DH10B_10reads.sam

    xenomappability --mapped_test_data tests/data/test_from_EcoliK12DH10B_10reads.sam > tests/data/test_from_EcoliK12DH10B_10reads.wig
    xenomappability --single_end_wiggle tests/data/test_from_EcoliK12DH10B_10reads.wig --sam_for_sizes tests/data/paired_end_testdata_human.sam`
	
Contributing to Xenomapper
==========================
Xenomapper is licensed under the GPLv3.  You are free to fork this repository under the terms of that license.  If you have suggested changes please start by raising an issue in the issue tracker.  Pull requests are welcome and will be included at the discretion of the author, but must have 100% test coverage.
Bug reports should be made to the issue tracker.  Difficulty in understanding how to use the software is a documentation bug, and should also be raised on the issue tracker with the tag `question` so your question and my response are easily found by others.


Citing Xenomapper
=================

Xenomapper is published in the Journal of Open Source Software.  Please cite the paper in academic publications [DOI:10.21105.joss.00018](http://dx.doi.org/10.21105.joss.00018).  Each release also has a Zenodo DOI identifier for each release.  In an ideal world this is what you would cite to indicate the code you use, and make everything more reproducible but academic credit is better served at the moment by the paper. Try and include the Zenodo DOI or a version number in your methods.  The DOI for the current release is [![DOI](https://zenodo.org/badge/11450/genomematt/xenomapper.svg)](https://zenodo.org/badge/latestdoi/11450/genomematt/xenomapper)

References
=================
Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359. http://bowtie-bio.sourceforge.net/bowtie2/

Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015 12:357-60. https://github.com/infphilo/hisat
