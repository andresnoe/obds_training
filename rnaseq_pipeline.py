'''
Dock string:
RNAseq pipeline

Start at fastq files
End with count matrix
'''

import gzip
from ruffus import *
from cgatcore import pipeline as P
import sys

@follows(mkdir("fastqc")) # makes folder in the directory in which you are running the function


@transform('*.fastq.gz', regex(r'(.*).fastq.gz'),r'fastqc/\1_fastqc.html')
def fastqc(infile, outfile):
    statement = "fastqc --nogroup -o fastqc %(infile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv="obds-py3")

@merge(fastqc, 'fastqc/multiqc_report.html')
def multiqc(infiles, outfile):
    statement = "multiqc -f -n %(outfile)s fastqc"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv="obds-py3")

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
