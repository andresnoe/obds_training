'''
run cellranger on 10X fastq files

'''
import re
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P
import sys

params = P.get_parameters("pipeline_cellranger.yml")  # define params here
samples = pd.read_csv("cellranger_samples.csv")
samples.set_index('name', inplace=True)
print(samples)

@follows(mkdir("count"))
@transform('data/*/.sample', regex(r'data/(.+)/.sample'), r'count/\1/outs/filtered_feature_bc_matrix.h5') #regex function from ruffus
def cellranger_count(infile, outfile):
    # python module looking for regular expression, group(1) is equiv to '\1'
    sampleid = re.search('data/(.+)/.sample', infile).group(1)
    print(sampleid)
    fastqs = samples['fastqs'][sampleid]
    cellnumber = samples['cells'][sampleid]
    chemistry = samples['chemistry'][sampleid]

    statement = '''cellranger count
    --id=               %(sampleid)s
    --transcriptome=    %(cellrangercount_transcriptome)s
    --fastqs=           %(fastqs)s
    --expect-cells=     %(cellnumber)s
    --chemistry=        %(chemistry)s
    --localcores=       %(cellrangercount_threads)s
    --localmem=         %(cellrangercount_memory)s
    > %(sampleid)s_standardout.log
    2> %(sampleid)s_standarderror.log
    && mv %(sampleid)s count/
    ''' # command line string statement, all options need to sit within these quotes, which p.run will send to the cluster as a job
        # Once cellranger runs, move the samples to count/ folder
        # If cellranger see count
    totalmem=str(params["cellrangercount_memory"])+'G'
    P.run(statement, job_queue='all.q',
          job_threads=params["cellrangercount_threads"],
          job_total_memory=totalmem,
          job_condaenv='obds-py3')

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
