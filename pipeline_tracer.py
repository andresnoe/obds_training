'''
Run TraCeR on all cells within a flow cell

'''
import re
import pandas as pd
from ruffus import *
from cgatcore import pipeline as P
import sys
import glob
import os

params = P.get_parameters("pipeline_tracer.yml")  # define params here

@transform("*.fastq.gz", regex(r'(.*).fastq.gz'), r'fastqc_results/\1_fastqc.html')
def fastqc(infile, outfile):
    ''' run fastqc on all files'''

    statement = '''fastqc --nogroup -o fastqc %(infile)s'''

    P.run(statement, job_queue=PARAMS["q"], job_threads=1, job_memory='2G', job_condaenv='tracer')

@follows(mkdir("multiqc_results"))
@merge(fastqc, "multiqc_results/fastqc_report.html")
def multiqc(infiles, outfile):
    ''' run multiqc to collect fastq stats '''

    statement = """export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastqc/ -f -n %(outfile)s"""
    P.run(statement, job_queue=PARAMS["q"])

# TRACER assemble FUNCTION HERE:
# What I've run previously that worked:
# nohup tracer assemble --ncores 16 --config_file /ifs/obds-training/sep20/andres/tracer/tracer/tracer.conf
# --loci A B --resource_dir /ifs/obds-training/sep20/andres/tracer/tracer/resources --species Hsap
# --quant_method kallisto --small_index --single_end --fragment_length 74 --fragment_sd 1
# P2_A4_S8_R1_001.fastq.gz 201119_test_cell 201119_test &

@tranform("*.fastq.gz", regex(r'(.*)_R[12]_001.fastq.gz'), r'%(assemble_output.directory)s/\1/filtered_TCR_seqs/filtered_TCRs.txt') #regex function from ruffus
def tracer_assemble(infiles, outfiles):

# /filtered_TCR_seqs/filtered_TCRs.txt
    '''run tracer assemble on fastq fastq.gz files
    Multiple outputs:
    For each cell, an /<output_directory>/<cell_name> directory will be created.
    This will contain 6 subdirectories:
    1) aligned_reads - Bowtie2 output
    2) Trinity_output - fasta files for each locus where contigs assembled
    3) IgBLAST_output -
    4) unfiltered_TCR_seqs - TCR sequences that were assembled prior to filtering
    5) expression_quantification - Kallisto/Salmon expression quantification
    6) filtered_TCR_seqs - recombinants filtered so 2 most expressed from each locus retained
    '''

    statement = '''rm -r %(outdir)s &&
    tracer assemble
    --ncores %(assemble_ncores)s
    --config_file %(tracer_options_config.file)s
    --loci %(tracer_options_loci)s
    --resource_dir %(tracer_options_resource.dir)s
    --species %(tracer_optionss_species)s
    --quant_method %(assemble_quant.method)s

    %(assemble_single.end)s
    --fragment_length %(assemble_fragment.length)s
    --fragment_sd %(assemble_fragment.sd)s
    %(assemble_small.index)s

    %(infiles_fastq1)s                     # These are wrong - look below
    OPTIONAL: %(infiles_fastq2)s           # These are wrong
    <cell_name>
    %(assemble_output.directory)s

     > %(sampleid)s_standardout.log
     2> %(sampleid)s_standarderror.log'''

    P.run(statement, job_queue=PARAMS["q"],
          job_threads=PARAMS["assemble_ncores"],
          job_memory=PARAMS["assemble_job.memory"],
          job_condaenv='tracer')



# nohup tracer summarise --config_file /ifs/obds-training/sep20/andres/tracer/tracer/tracer.conf
# --loci A B --resource_dir /ifs/obds-training/sep20/andres/tracer/tracer/resources
# --species Hsap --keep_invariant --no_networks
# /ifs/obds-training/sep20/andres/VAC066_sc_bcl/190621_NB501183_0664_AHT7MHBGXB/Data/Intensities/BaseCalls/$/test_output &


@split("*.fastq.gz", regex(r'(.*)_[12].fastq.gz'), r'tracer_results/%(batch)s\1') #regex function from ruffus
# WANT TO APPEND 'BATCH 1' to BATCH 1 cell names - can I parameterise this?
def tracer_summarise(infiles, outfiles):

    '''run tracer summarise on output of tracer_assemble
    Multiple outputs to <input_dir>/filtered_TCR_summary
    1) TCR_summary.txt
    2) recombinants.txt
    3) reconstructed_lengths_TCR[A|B]
    4) clonotype_sizes
    5) and 6) only enabled in Graphviz installed
    (see https://github.com/Teichlab/tracer)
    '''

    statement = '''
    tracer summarise
    --config_file %(tracer_options_config.file)s
    --loci %(tracer_options_loci)s
    --resource_dir %(tracer_options_resource.dir)s
    --species %(tracer_optionss_species)s

--keep_invariant --no_networks
    <input_dir>
    '''
    P.run(statement, job_queue=PARAMS["q"],
          job_threads= 6,
          job_memory= 2G,
          job_condaenv='tracer')



if __name__ == "__main__":
    sys.exit(P.main(sys.argv))



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
