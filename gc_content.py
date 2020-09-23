#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:51:06 2020

@author: andresnoe

Calculate GC content
"""



# ENST00000530893.2 BRCA2-003 cdna:protein_coding

sequence = '''GGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCGCCTCGGGTGTCTTTTGC
GGCGGTGGGTCGCCGCCGGGAGAAGCGTGAGGGGACAGATTTGTGACCGGCGCGGTTTTT
GTCAGCTTACTCCGGCCAAAAAAGAACTGCACCTCTGGAGCGGACTTATTTACCAAGCAT
TGGAGGAATATCGTAGGTAAAAATGCCTATTGGATCCAAAGAGAGGCCAACATTTTTTGA
AATTTTTAAGACACGCTGCAACAAAATTTAGGACCAATAAGTCTTAATTGGTTTGAAGAA
CTTTCTTCAGAAGCTCCACCCTATAATTCTGAACCTGCAGAAGAATCTGAACATAAAAAC
AACAATTACGAACCAAACCTATTTAAAACTCCACAAAGGAAACCATCTTATAATCAGCTG
GCTTCAACTCCAATAATATTCAAAGAGCAAGGGCTGACTCTGCCGCTGTACCAATCTCCT
GTAAAAGAATTAGATAAATTCAAATTAGACTTAGGAAGGAATGTTCCCAATAGTAGACAT
AAAAGTCTTCGCACAGTGAAAACTAAAATGGATCAAGCAGATGATGTTTCCTGTCCACTT
CTAAATTCTTGTCTTAGTGAAAGTCCTGTTGTTCTACAATGTACACATGTAACACCACAA
AGAGATAAGTCAGTGGTATGTGGGAGTTTGTTTCATACACCAAAGTTTGTGAAGGGTCGT
CAGACACCAAAACATATTTCTGAAAGTCTAGGAGCTGAGGTGGATCCTGATATGTCTTGG
TCAAGTTCTTTAGCTACACCACCCACCCTTAGTTCTACTGTGCTCATAGTCAGAAATGAA
GAAGCATCTGAAACTGTATTTCCTCATGATACTACTGCTAATGTGAAAAGCTATTTTTCC
AATCATGATGAAAGTCTGAAGAAAAATGATAGATTTATCGCTTCTGTGACAGACAGTGAA
AACACAAATCAAAGAGAAGCTGCAAGTCATGGATTTGGAAAAACATCAGGGAATTCATTT
AAAGTAAATAGCTGCAAAGACCACATTGGAAAGTCAATGCCAAATGTCCTAGAAGATGAA
GTATATGAAACAGTTGTAGATACCTCTGAAGAAGATAGTTTTTCATTATGTTTTTCTAAA
TGTAGAACAAAAAATCTACAAAAAGTAAGAACTAGCAAGACTAGGAAAAAAATTTTCCAT
GAAGCAAACGCTGATGAATGTGAAAAATCTAAAAACCAAGTGAAAGAAAAATACTCATTT
GTATCTGAAGTGGAACCAAATGATACTGATCCATTAGATTCAAATGTAGCAAATCAGAAG
CCCTTTGAGAGTGGAAGTGACAAAATCTCCAAGGAAGTTGTACCGTCTTTGGCCTGTGAA
TGGTCTCAACTAACCCTTTCAGGTCTAAATGGAGCCCAGATGGAGAAAATACCCCTATTG
CATATTTCTTCATGTGACCAAAATATTTCAGAAAAAGACCTATTAGACACAGAGAACAAA
AGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGCCTACCAAAA
TCAGAGAAGCCATTAAATGAGGAAACAGTGGTAAATAAGAGAGATGAAGAGCAGCATCTT
GAATCTCATACAGACTGCATTCTTGCAGTAAAGCAGGCAATATCTGGAACTTCTCCAGTG
GCTTCTTCATTTCAGGGTATCAAAAAGTCTATATTCAGAATAAGAGAATCACCTAAAGAG
ACTTTCAATGCAAGTTTTTCAGGTCATATGACTGATCCAAACTTTAAAAAAGAAACTGAA
GCCTCTGAAAGTGGACTGGAAATACATACTGTTTGCTCACAGAAGGAGGACTCCTTATGT
CCAAATTTAATTGATAATGGAAGCTGGCCAGCCACCACCACACAGAATTCTGTAGCTTTG
AAGAATGCAGGTTTAATATCCACTTTGAAAAAGAAAACAAATAAGTTTATTTATGCTATA
CATGATGAAACATCTTATAAAGGAAAAAAAA
'''

guanine = 0
cytosine = 0
adenine = 0
thymine = 0
unknown = 0

for base in sequence:
    if base == "G":
        guanine += 1
    elif base == "C":
        cytosine += 1
    elif base == "T":
        thymine += 1
    elif base == "A":
        adenine += 1
    elif base == '\n':
        pass
    else:
        unknown += 1

gc_content = ((guanine+cytosine)/(guanine+cytosine+adenine+thymine+unknown)*100)

print(f'GC content is {gc_content}')