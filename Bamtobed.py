# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:18:45 2020

Bamtobed pysam
# BED file format 
1) chromosome 
2) chrom start (0 based)) 
3) chrom end 
4) Name 
5) Score 
6) Strand 
@author: surface
"""
# Convert Sam to BED using Pysam
# read in Sam file in Pysam
# iterate over the sam file line by line
# Create the bed file
# import pysam


import pysam
samfilepath = ""
samfile = pysam.AligmentFile(samfilepath,"r")
samfile_iter = samfile.fetch()

for alignment in samfile_iter:
    if alignment.is_paired:
        print(f"{alignment.reference_name}\t{alignment.query_alignment_start}\t{alignment.query_alignment_end}\t{alignment.query_name}\t{alignment.mapping_quality}\t.\n")
        print(alignment)
        


