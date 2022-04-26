import argparse
import subprocess
import time
from helpers import check_fail

def run_diamond(index, in_path, out_path, threads, outfmt):
    diamond_command = "diamond" + " blastx --db " + index +\
                        " --query " + in_path + " --sensitive --max-target-seqs 1 --outfmt " + str(outfmt) +\
                        " --masking 0 -c 1 -s 1 -b 0.75" +\
                        " --threads " + str(threads) +\
                        " --out " + out_path
    new_command = subprocess.run(diamond_command, shell=True)
    if check_fail("diamond", new_command, []) is True: return False

def diamond(args: argparse.Namespace):
    dirpath = args.dirpath

    new_fwd = dirpath + "/new_fwd.fq"
    new_rev = dirpath + "/new_rev.fq"

    new_contigs = contig_path + "/final_contigs.fq"
    contig_path = dirpath + "/megahit_out"

    # now lets align reads
    # need to merge paired end reads first though
    merged_pe = dirpath + "/merged_reads.fq"
    merge_command = "seqtk mergepe " + new_fwd + " " + new_rev + " > " + merged_pe
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("seqtk mergepe", new_command, []) is True: return False 

    diamonc_combined_file = dirpath + "/diamond_combined_file.fq"
    subprocess.run("cat " + new_contigs + " > " + diamonc_combined_file, shell=True)
    subprocess.run("cat " + merged_pe + " >> " + diamonc_combined_file, shell=True)

    nr_combined_file = dirpath + "/nr_combined_file"
    start = time.time()
    if run_diamond(args.diamond_index, diamonc_combined_file, nr_combined_file, args.threads, 6) == False: return None
    end = time.time()
    print("contig alignment against nr took: " + str(end - start))