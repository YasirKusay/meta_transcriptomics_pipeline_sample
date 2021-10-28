import argparse
import subprocess
import tempfile
import shutil
import os
#from main_pipeline.helpers import check_command_exists, check_fail, generate_temp_file
from helpers import check_command_exists, check_fail, generate_temp_file


def run_pipeline(args: argparse.Namespace):

    # setting up a temporary dir for all intermediate files generated
    #dirpath = tempfile.mkdtemp()

    # lets instead check if pipeline works ok :)
    dirpath = args.dirpath

    ##################### FASTP ########################

    path_command = check_command_exists('fastp')

    # generating temp files, to be removed later
    #out1 = "/srv/scratch/z5215055/HONS/res/out1.fastq"
    #out2 = "/srv/scratch/z5215055/HONS/res/out2.fastq"
    generated_files = []

    out1 = generate_temp_file("fastq", dirpath)
    out2 = generate_temp_file("fastq", dirpath)
    generate_temp_file.append(out1, out2)

    fastp_path = "fastp"
    fastp_command = fastp_path +\
                    " --in1 " + args.inp1 +\
                    " --in2 " + args.inp2 +\
                    " --out1 " + out1 +\
                    " --out2 " + out2 +\
                    " --qualified_quality_phred  " + args.qualified_quality_phred +\
                    " --unqualified_percent_limit " + args.unqualified_percent_limit +\
                    " --length_required " + args.length_required +\
                    " --thread " + str(args.threads)
                    # need to consider adapters, should we give the user a chance to add adatpers?

    new_command = subprocess.run(fastp_command, shell=True, capture_output=True)
    if check_fail(new_command, [out1, out2]) is False: return None

    print("PART 1 DONE FASTP")
    print(new_command.returncode)
    print('out: ', new_command.stdout.decode())
    print('err: ', new_command.stderr.decode()) 

    #################################### SORTMERNA ############################
    aligned = generate_temp_file("", dirpath)
    other = generate_temp_file("", dirpath)
    
    path_command = check_command_exists('sortmerna')
    sortmerna_path = 'sortmerna'
    sortmerna_command = sortmerna_path +\
                    " --ref  " + args.sortmerna_index +\
                    " --aligned " + aligned +\
                    " --other " + other +\
                    " --fastx " +\
                    " --reads " + out1 + " --reads " + out2 +\
                    " -threads " + str(args.threads) +\
		            " --out2 TRUE " +\
		            " --paired_out TRUE" +\
                    " --best 1 " # 1 = all high candidate reference sequences will be searched for alignments

    new_command = subprocess.run(sortmerna_command, shell=True, capture_output=True)
    if check_fail(new_command, [out1, out2]) is False: return None
    os.remove(aligned + "_fwd.fastq", aligned + ".log", aligned + "_rev.fastq")
    generated_files.append(other + "_fwd.fastq", other + "_rev.fastq")

    #if check_fail(new_command, [file1, file2]) is False: return None
    print("PART 2 DONE SORTME")
    print(new_command.returncode)
    print('out: ', new_command.stdout.decode())
    print('err: ', new_command.stderr.decode())


    ############################# SNAP HUMAN #############################
    path_command = check_command_exists('snap-aligner')
    human_out = generate_temp_file("sam", dirpath) # stores human output in sam file, can be explored further to align reads, etc.
    snap_path = 'snap-aligner'
    snap_human_command = snap_path + " paired " + args.snap_human_index + other + "_fwd.fastq " + other + "_rev.fastq " +\
                    " -o " + human_out + " -t " + str(args.threads)
    new_command = subprocess.run(snap_human_command, shell=True, capture_output=True)
    if check_fail(new_command, [human_out]) is False: return None

    path_command = check_command_exists('samtools')
    human_subtract_1 = generate_temp_file("fastq", dirpath)
    human_subtract_2 = generate_temp_file("fastq", dirpath)
    human_spare = generate_temp_file("fastq", dirpath)

    # retrieving only human reads
    samtools_path = "samtools"
    samtools_human_subtrack_command = samtools_path + " fastq  -f 4 -@ " + str(args.threads) +\
                        " -1 " + human_subtract_1 +\
                        " -2 " + human_subtract_2 +\
                        " -s " + human_spare +\
                        human_out

    new_command = subprocess.run(snap_human_command, shell=True, capture_output=True)
    if check_fail(new_command, [human_subtract_1 + ".fastq", human_subtract_2 + ".fastq", human_spare + ".fastq"]) is False: return None
    generated_files.append(human_subtract_1 + ".fastq", human_subtract_2 + ".fastq", human_spare + ".fastq")

    #################### MEGAHIT ###########################
    megahit_path = "megahit"
    contig_path = dirpath + "_out"
    megahit_command = megahit_path + " -1 " + human_subtract_1 +\
                        " -2 " + human_subtract_2 +\
                        " -o " + contig_path # is an output directory 

    contigs = contig_path + "/final_contigs.fa"

    # we can now align the contigs to the databases
    # lets do them sequentially for now
    snap_contigs = generate_temp_file("sam", dirpath)
    snap_contig_command = snap_path + " single " + contigs +\
                    " -o " + snap_contigs + " -t " + str(args.threads)
    new_command = subprocess.run(snap_contig_command, shell=True, capture_output=True)
    if check_fail(new_command, [human_out]) is False: return None

    # now lets align reads
    diamond_path = "diamond"
    diamond_contigs = generate_temp_file("sam", dirpath)
    diamond_command = diamond_path + " blastx -db " + args.diamond_index +\
                        " --query " + contigs + " --sensitive --max-target-seqs 1 --outfmt 101"\
                        " --threads " + args.threads +\
                        " --out " + diamond_contigs

    # we are done, lets remove the temp directory
    #shutil.rmtree(dirpath)