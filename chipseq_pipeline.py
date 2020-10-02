#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChIP-seq peakcalling pipeline

input:
    - chipseq bam files & their associated input files

outputs:
    - filtered bamfiles
        - remove non uniquely mapped reads
        - remove duplicates
        - remove reads in encode blacklist
    - multiqc html report
    - peakfiles
        - narrowPeaks (from macs2)
        - summits.bed (from macs2)
        - merged replicates
    - bigwigs
        - cpm normalised
cd
"""

import sys
from ruffus import *
from cgatcore import pipeline as P

PARAMS = P.get_parameters("chipseq_peakcalling_pipeline.yml")


# Filter bamfile ​remove unmapped read & MAPQ quality 20 ​
@follows(mkdir('filtered_bam'))
@transform("*.bam", regex(r'(.*).bam'), r'filtered_bam/\1.bam')
def filter_bam(infile, outfile):
    statement = """
    samtools view
    -b -h
    %(filterbam_opts)s
    %(infile)s
    > %(outfile)s
    && samtools index %(outfile)s
    """
    P.run(statement, job_queue=PARAMS["q"],
          job_threads=PARAMS["threads"],
          job_memory=PARAMS["job_memory"],
          job_condaenv='obds-py3')

# remove duplicate reads ​​
@transform(filter_bam, regex(r'filtered_bam/(.*).bam'), r'filtered_bam/\1.dedup.bam')
def filtered_duplicate(infile, outfile):
    statement="""
       picard MarkDuplicates
       -Xmx%(picard_memory)sg
       I=%(infile)s
       O=%(outfile)s
       M=%(outfile)s.metrics
       REMOVE_DUPLICATES=%(picard_remove_dups)s"""
    final_memory = str(int(PARAMS['picard_memory']) + 2) + 'g'
    P.run(statement, job_queue=PARAMS["q"],
                     job_memory=final_memory,
                     job_threads=PARAMS["threads"],
                     job_condaenv='obds-py3')

# Remove reads in blacklisted regions ​
@transform(filtered_duplicate, regex(r'filtered_bam/(.*).dedup.bam'), r'filtered_bam/\1.dedup.blacklist.bam')
def filtered_blacklist(infile, outfile):
    statement= """
        alignmentSieve
        -bl %(blacklist)s
        -b %(infile)s
        -p %(deeptools_threads)s
        -o %(outfile)s"""
    P.run(statement, job_queue=PARAMS["q"],
                    job_threads=PARAMS["deeptools_threads"],
                    job_condaenv='peaktools_env')
                    #change env for Deeptools as this isn't in obds-py3 env.


# Call peaks (MACS2)​
#https://regex101.com/
@follows(mkdir('peaks'))
@transform(filtered_blacklist, regex(r'))
def macs2(infile, outfile):
statement="""macs2 callpeak
            -t %(chipfile)s
            -c %(inputfile)s
            -n %(outfile)s
            --outdir %(peaks)s
            -f %(BAM)s
            """)

# Compare peaklist across replicates (Bedtools)

# Generate Bigwigs (Deeptools)

# Merge GR peaks to get consensus peak list (Bedtools)




if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
