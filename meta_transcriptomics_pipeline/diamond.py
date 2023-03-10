import argparse
import subprocess
import time
from meta_transcriptomics_pipeline.helpers import check_fail

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

    nr_combined_file = dirpath + "/nr_combined_file"
    combined_file = dirpath + "/combined_file.fq"

    start = time.time()
    if run_diamond(args.diamond_index, combined_file, nr_combined_file, args.threads, 6) == False: return None
    end = time.time()
    print("contig alignment against nr took: " + str(end - start))