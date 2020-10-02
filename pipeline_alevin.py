'''
run salmon alevin on 10X fastq files

'''
import re
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P
import sys
import glob
import os

params = P.get_parameters("pipeline_alevin.yml")  # define params here
samples = pd.read_csv("cellranger_samples.csv")
samples.set_index('name', inplace=True)
print(samples)

def get_gex_fastq(dir):
    '''Docstring'''
    fastq1_pattern = params["pattern"]["fastq1"]
    fastq1_glob = f"{dir}/*{fastq1_pattern}*"
    fastq1 = glob.glob(fastq1_glob)
    if len(fastq1) == 0:
        raise OSError(f"No file matched pattern: {fastq1_glob}")
    fastq2 = [file.replace(params["pattern"]["fastq1"], params["pattern"]["fastq2"]) for file in fastq1]
    #For each file of fastq1, replace the name of fastq2 to be the same
    for file in fastq2:
        if not os.path.exists(file):
            raise OSError(f"Paired file not found: {file}")
    return {'fastq1' : fastq1, 'fastq2' : fastq2 } # returns a dictionary

@follows(mkdir("alevin"))
@transform('data/*/.sample', regex(r'data/(.+)/.sample'), r'alevin/\1/alevin/quants_mat.gz') #regex function from ruffus
def salmon_alevin(infile, outfile):
    # python module looking for regular expression, group(1) is equiv to '\1'
    outdir = outfile.replace("/alevin/quants_mat.gz", "")
    sampleid = re.search('data/(.+)/.sample', infile).group(1)
    print(sampleid)
    fastqs = samples['fastqs'][sampleid]
    cellnumber = samples['cells'][sampleid]
    chemistry = samples['chemistry'][sampleid]

    fastq_dict = get_gex_fastq(fastqs) # calls function above to find fastq files
    print(fastq_dict)

    infiles_fastq1 = " ".join(fastq_dict["fastq1"])
    infiles_fastq2 = " ".join(fastq_dict["fastq2"])

    chemistry = samples['chemistry'][sampleid]
    if chemistry == "SC3Pv2":
        chemistry_option = "--chromium"
    elif chemistry == "SC3Pv3":
        chemistry_option = "--chromiumV3"
    else:
        raise NameError('Invalid chemistry.')

    statement = '''rm -r %(outdir)s &&
    salmon alevin
     -l ISR
     -1 %(infiles_fastq1)s
     -2 %(infiles_fastq2)s
     %(chemistry_option)s
     -i %(index_file)s
     -p %(threads)s
     -o %(outfile)s
     --tgMap %(tgMap_file)s
     --expectCells %(cellnumber)s
     > %(sampleid)s_standardout.log
     2> %(sampleid)s_standarderror.log'''

    P.run(statement, job_queue='all.q',
          job_threads=params["threads"],
          job_memory=params["job_memory"],
          job_condaenv='obds-py3')

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
