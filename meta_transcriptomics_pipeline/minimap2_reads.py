import argparse
import subprocess
import time
import os
from meta_transcriptomics_pipeline.helpers import check_fail
from meta_transcriptomics_pipeline.paf2blast6 import paf2blast6

def minimap2_reads(args: argparse.Namespace):
    start = time.time()
    dirpath = args.dirpath
    bigger_1 = dirpath + "/bigger_1.fq"
    bigger_2 = dirpath + "/bigger_2.fq"
    minimap2_path = "minimap2"
    minimap2_long_reads_out = dirpath + "/minimap2_lr_out.paf"
    tempdir = dirpath + "/tempdir2"
    os.mkdir(tempdir)
    minimap2_command = minimap2_path + " -c -x sr -t " + str(args.threads) + " --split-prefix " + tempdir + "/temp" + " " + args.minimap2_index + " " + bigger_1 + " " + bigger_2 + " > " + minimap2_long_reads_out
    new_command = subprocess.run(minimap2_command, shell=True)
    if check_fail(minimap2_path, new_command, []) is True: return None

    end = time.time()
    print("minimap2 alignment 2 against nt took: " + str(end - start))

    paf2blast6(minimap2_long_reads_out, dirpath)
    paf2blast_out = dirpath + "/minimap2_lr_out_frompaf.m8"