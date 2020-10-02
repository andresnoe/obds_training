''' Script to split-transform-merge as pipeline


Split file by chrom
Count the number of transcripts per chrom
Read all count files
    Calculate the average
    Write to files

    '''

from ruffus import *
from cgatcore import pipeline as P
import gzip



@split(‘test.gtf’, ‘chr*.gtf’)
def split_chrom(infile, outfiles):
    with gzip.open(infile, 'r') as inf:
        last_chrom=""
        for line in inf:
            chrom = line.split()[0]
            if chrom == last_chrom:
                outfile.write(line)
            else:
                if last_chrom != "":
                    outf.close
                outfile = chrom + ".gtf.gz"
                outf = gzip.open(outfile, 'wt')
                outf.write(line)
                last_chrom=chrom
                print(chrom)
# split_chrom("test.gtf.gz", "")

@transform('chr*.gtf.gz', suffix('.gtf.gz'),'.count')
def count_genes(infile, outfile):
    statement = "wc -l %(infile)s > %(outfile)s"
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv="obds-py3")


if __name__ == "__rain__":
    sys.exit(P.main(sys.argv))
