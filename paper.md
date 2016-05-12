---
title: 'Xenomapper: Mapping reads in a mixed species context'
tags:
  - bioinformatics
  - short read mapping
  - multispecies
  - xenograft
authors:
 - name: Matthew J. Wakefield
   orcid: 0000-0001-6624-4698
   affiliation: The Walter and Eliza Hall Institute
date: 12 May 2016
bibliography: paper.bib
---

# Summary

Xenomapper is a utility for post processing mapped DNA sequencing reads that have been aligned to a primary genome and a secondary genome, and binning reads into species specific, multimapping in each species, unmapped and unassigned categories. It can be used on single end or paired end sequencing data across a wide range of genomics methods including RNAseq. In paired end data evidence of sequence specificity for either read will be used to assign both reads.

Use cases include xenografts of human cancers and host pathogen interactions.

Xenomapper is most effective with mapped reads that include an XS or ZS score that gives the mapping score of the next best read. These include Bowtie2 [@BOWTIE2] and HISAT [@HISAT].

This work builds upon similar a approaches by Rossello et. al. [@ROSSELLO] with a more rigorous implementation and extensions for paired end data and exon aware aligners.

# References

