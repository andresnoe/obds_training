#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 10:34:33 2020

overlap using chromStart chromend

use brain bedfile as reference and comparing gut bedfile to every chromStart and chromEnd

First attempt: 
    Comparing B (GUT) start to A (BRAIN) chrstart and end
    Comparing B end to A chrstart and end
    
TO DO:
1) Print number of intersects

2) Print overlapping coordinates from one file
    
3) Include arg.parse() etc 

4) zip at the end

5) add log files l...


@author: Jennifer, Hebe, Sophie and Andrés
"""




# a = /Users/jennifer/Desktop/week2/week2/test_overlap/brain_dnase1_chr21.bed
# b = /Users/jennifer/Desktop/week2/week2/test_overlap/gut_dnase1_chr21.bed

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input1', '-i1', dest='brainbedfilepath', help='input path file')
parser.add_argument('--input2', '-i2', dest='gutbedfilepath', help='input path file')
parser.add_argument('--output', '-o', dest='outputfilepath', help='output count path file')
args = parser.parse_args()
# parser.add_argument('--padding', '-p', dest='bedfilepad', type= int, default=0, help='number of bases added')
# args = parser.parse_args()


# with open(args.samfilepath, 'r') as samfile:
    # with gzip.open(args.bedfilepath, 'wt') as bedfile:

overlap_count=0   

with open(args.brainbedfilepath,'r') as brainbedfile, open(args.gutbedfilepath, 'r') as gutbedfile:
    with open(args.outputfilepath, 'wt') as outputfile:
        for line in brainbedfile:
            col_brain = line.split()
            chromStart_brain = int(col_brain[1])
            chromEnd_brain = int(col_brain[2])
            print(chromStart_brain,chromEnd_brain)
            for line in gutbedfile:
                col_gut = line.split()
                chromStart_gut = int(col_gut[1])
                chromEnd_gut = int(col_gut[2])
                print(chromStart_gut,chromEnd_gut)
                if (chromStart_brain >= chromStart_gut) and (chromStart_brain <= chromEnd_gut):
                    overlap_count += 1
                    # print(overlap_count)
                    outputfile.write(f'{overlap_count}\t{chromStart_brain}\t{chromEnd_brain}\n')

                
                
            # print(chromStart_gut) 
            # print(chromStart_brain)
            # print(chromEnd_brain)
                  
 