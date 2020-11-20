'''
Input = fastq  files

Steps required
- 1_fastqc
- 2_multiqc

To run full pipeline: python pipeline_multiqc.py make
To check up to and including a point:
python pipeline_multiqc.py show multiqc -v 5

'''

from ruffus import *
from cgatcore import pipeline as P
import sys

PARAMS = P.get_parameters("kallisto_pipeline_a.yml")

@follows(mkdir("fastqc"))
@transform("*.fastq.gz", regex(r'(.*).fastq.gz'), r'fastqc/\1_fastqc.html')
def fastqc(infile, outfile):
    ''' run fastqc on all files'''

    statement = '''fastqc --nogroup -o fastqc %(infile)s'''
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')

@follows(mkdir("reports"))
@merge(fastqc, "reports/fastqc_report.html")
def multiqc(infiles, outfile):
    ''' run multiqc to collect fastq stats '''

    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastqc/ -f -n %(outfile)s"""
    P.run(statement, job_queue=PARAMS["q"])

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
