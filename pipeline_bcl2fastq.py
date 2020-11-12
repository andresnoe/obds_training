'''
Simple pipeline to use Illumina's bcl2fastq cli
https://bioinformatics.cvr.ac.uk/how-to-demultiplex-illumina-data-and-generate-fastq-files-using-bcl2fastq/
https://www.outils.genomique.biologie.ens.fr/samplesheetvalidator/
https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf

'''
# NEED TO: change samplesheet on server
# Trying this command: nohup bcl2fastq --runfolder-dir . -p 16 --output-dir /190621_NB501183_0664_AHT7MHBGXB/fastq_files &
# Got this error: 2020-11-09 20:53:54 [2b1c84f446a0] ERROR: bcl2fastq::common::Exception: 2020-Nov-09 20:53:54: Permission denied (13): /home/shahlab/miniconda2/conda-bld/bcl2fastq_1497316272077/work/bcl2fastq/src/cxx/lib/common/FileSystem.cpp(40): Throw in function void bcl2fastq::common::createDirectories(std::vector<boost::filesystem::path>)
#Dynamic exception type: boost::exception_detail::clone_impl<bcl2fastq::common::IoError>
# std::exception::what: Failed to create directory /190621_NB501183_0664_AHT7MHBGXB/fastq_files/Stats

# TRied this
# nohup bcl2fastq --runfolder-dir . -p 16 &
# Worked  but split across  lanes.
# Specify  --no-lane-splitting next time

#Tried this:
# nohup bcl2fastq --runfolder-dir . -p 32 -r 4 -w 4 --no-lane-splitting next time &


from ruffus import *
from cgatcore import pipeline as P
import sys

PARAMS = P.get_parameters("pipeline_bcl2fastq.yml")

@follows(mkdir("fastq_files"))
@transform("*.fastq.gz", regex(r'(.*).fastq.gz'), r'fastqc/\1_fastqc.html') # NOT YET RIGHT
def bcl2fastq(infile, outfile):
    ''' Run bcl2fastq on one flow cell'''

    statement = '''
                    bcl2fastq
                    --runfolder-dir `pwd`
                    --output-dir %(outfile)s
                    --sample-sheet SampleSheet.csv
                    --use-bases-mask "Y*,I8,I8"
                    --barcode-mismatches 1
                    --ignore-missing-bcl
                    --ignore-missing-control
                    --loading-threads
                    --processing-threads
                    --writing-threads'''
    P.run(statement, job_queue='all.q', job_threads=1, job_memory='2G', job_condaenv='obds-py3')



if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


#    @posttask(touch_file(os.path.join(runs_scratch_dir,'fastqs','completed')))
#def bcl2fastq_conversion(run_directory, completed_flag):
#    """ Run bcl2fastq conversion and create fastq files in the run directory"""
#    out_dir = os.path.join(runs_scratch_dir,'fastqs')
#    interop_dir = os.path.join(out_dir,'InterOp')
#
#    # r, w, d, and p specify numbers of threads to be used for each of the concurrent subtasks of the conversion (see bcl2fastq manual)
#    args = "-R {indir} -o {outdir} --interop-dir={interopdir} -r1 -w1 -d2 -p4 \
#            ".format(indir=run_directory, outdir=out_dir, interopdir=interop_dir)
#    if options.run_on_bcl_tile != None:
#        args += " --tiles %s" % options.run_on_bcl_tile
