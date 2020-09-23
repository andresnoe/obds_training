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
import gzip

# path for samfile
# samfilepath = "/Users/andresnoe/obds_sep20/working_directory/ERR1755082.test.sam"
# bedfilepath = "/Users/andresnoe/obds_sep20/working_directory/output.bed"

parser = argparse.ArgumentParser(description=("Script to convert SAM to BED format"))
parser.add_argument('--input' , '-i', dest='samfilepath',
                    help='input file path')
parser.add_argument('--output','-o', dest='bedfilepath',
                    help='output file path')
parser.add_argument('--pad','-p', dest='bedfilepad',
                    type =int,
                    default = 0,
                    help='number of bases added to start and end of bed file interval')
parser.add_argument('--fragment', '-f',less dest='fragment',
                    action = 'store_true',
                    default=False,
                    help="Boolean flag: the output bedfile provides the coordinates of the sequenced FRAGMENTS rather than reads (default = FALSE)")
args=parser.parse_args()

with open (args.samfilepath, 'r') as samfile:
    #Here samfile is the variable where the samfile will be held in python
    with gzip.open(args.bedfilepath, 'wt') as bedfile:
        for line in samfile:
            if line[0] == '@':
                pass
                
            else:
                # 
                col = line.split()
                # chromosome, (col 3 in SAM)
                chrom = col[2]
                if chrom == "*":
                    continue
                # start position in the chromosome, (col 4 but 0 based)
                # end position in the chromosome, (col 3 + len(col 10))
                if args.bedfilepad == 0 :
                    start = (int(col[3])-1)
                    end = int(col[4]) + int(len(col[9]))
                elif args.bedfilepad >0:
                    start = (int(col[3]) -1 -args.bedfilepad)
                    # Want to SUBTRACT THE PAD from the start site
                    end = int(col[4]) + int(len(col[9])) + args.bedfilepad
                    # Want to ADD THE PAD to the end site
                # name, (col 1)
                name = col[0]
                # score, (col 5)
                score = int(col[4])
                if score == 0:
                    continue
                # strand for a set of genomic positions (hard code ".")
                strand = '.'
                if args.fragment:
                # TLEN column in SAM file to avoid using stranded pairs twice
                    tlen = int(col[8])
                    if tlen>0:
                        end = (int(col[3]) + tlen)
                        bedfile.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n')
                else:
                    bedfile.write(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n')


