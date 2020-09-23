#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:16:12 2020

@author: andresnoe

Script designed to convert SAM file to BED file


SAM FILE:
    http://samtools.github.io/hts-specs/SAMv1.pdf
    page 6

BED FILE:
    A 6 column BED file contains:
        chromosome, (col 3 in SAM)
        start position in the chromosome, (col 4 but 0 based)
        end position in the chromosome, (col 3 + len(col 10))
        name, (col 1)
        score, (col 5)
        strand for a set of genomic positions (hard code ".")

"""
import argparse

# path for samfile
# samfilepath = "/Users/andresnoe/obds_sep20/working_directory/ERR1755082.test.sam"
# bedfilepath = "/Users/andresnoe/obds_sep20/working_directory/output.bed"

parser = argparse.ArgumentParser()
parser.add_argument('--input' , '-i', dest='samfilepath',
                    help='input file path')
parser.add_argument('--output','-o', dest='bedfilepath',
                    help='output file path')
args=parser.parse_args()

with open (args.samfilepath, 'r') as samfile:
    #Here samfile is the variable where the samfile will be held in python
    with open(args.bedfilepath, 'w') as bedfile:
        for line in samfile:
            if line[0] == '@':
                pass
                
            else:
                # 
                col = line.split()
                # chromosome, (col 3 in SAM)
                chrom = col[2]
                # start position in the chromosome, (col 4 but 0 based)
                start = (int(col[3])-1)
                #         end position in the chromosome, (col 3 + len(col 10))
                end = int(col[4]) + int(len(col[9]))
                #         name, (col 1)
                name = col[0]
                #         score, (col 5)
                score = col[4]
                #         strand for a set of genomic positions (hard code ".")
                strand = '.'
                bedfile.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n')

