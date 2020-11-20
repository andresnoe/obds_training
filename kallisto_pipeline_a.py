'''
Kallisto pipeline team A


Input = fastq  files

Steps required
- 1_fastqc
- 2_multiqc
- 3_kallisto_index
- 4_kallisto_quant
- 5_sleuth_differential expression (R)
- 6_pathway_analysis (R)

To run full pipeline: python kallisto_pipeline_a.py make
To check up to and including a point:
python kallisto_pipeline_a.py show kallisto_quant -v 5

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

@follows(fastqc, multiqc)
@follows(mkdir("kallisto_index_folder"))
# Input is *.fa.gz ---> is this the same as *.fasta.gz?
@transform(PARAMS["index_fasta"], regex(r'.*/(.*).fa.gz'), r'kallisto_index_folder/\1.idx')
def kallisto_index(infile, outfile):
    ''' generate index from genome '''
    statement = "kallisto index --index=%(outfile)s --kmer-size=31 %(infile)s"
    P.run(statement, job_queue=PARAMS["q"], job_threads=1, job_memory='8G', job_condaenv='obds-py3')

# ls -1 *test.fastq.gz
# ERR1755082_1_test.fastq.gz
# ERR1755082_2_test.fastq.gz
#
# ERR1755082_[12]_test.fastq.gz
#
# (.*)_[12]_test.fastq.gz

#regex(r'(.*)_[12](.*).fastq.gz')

# Kallisto manual: https://pachterlab.github.io/kallisto/manual
@follows(mkdir("kallisto_quant_folder"))
@follows(kallisto_index)
@collate("*.fastq.gz", regex(r'(.*)_[12].fastq.gz'), r'kallisto_quant_folder/\1')
# @collate("*.fastq.gz", regex(r'(.*)_[12](.*).fastq.gz'), r'bam/\1.bam')
def kallisto_quant(infiles, outfile):
    #read1 = '*1_test.fastq.gz'
    #read2 = '*2_test.fastq.gz'
    #fastq_infiles = infiles[0:-1]
    fastq_infiles = " ".join(infiles)
    '''run kallisto pseudo on fastq files
        3 outputs:  abundances.h5 HDF5 binary file containing run info
                    abundances.tsv of abundance estimates
                    run_info.json containing information about the run'''

    statement = '''
                kallisto quant
                --index=%(kallisto_quant_index)s
                --output-dir=%(outfile)s
                %(kallisto_options_bootstrap)s
                %(kallisto_options_strandness)s
                %(kallisto_options_seed)s
                %(fastq_infiles)s
                '''
    P.run(statement, job_queue=PARAMS["q"], job_threads=16, job_memory='2G', job_condaenv='obds-py3')


#Need to use stdout and stderror for next function
@follows(mkdir("kallisto_qc"))
@merge(kallisto_quant, r"kallisto_qc/fastqc_report.html")
def kallisto_qc(infiles, outfile):
    ''' run multiqc to collect fastq stats after kallisto_quant '''
    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc --force --filename %(outfile)s kallisto_quant_folder"""
    P.run(statement, job_queue=PARAMS["q"], job_threads=1, job_memory='2G', job_condaenv='obds-py3')


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
