import argparse
import subprocess
import time
import os
import shutil
from meta_transcriptomics_pipeline.helpers import check_fail
from meta_transcriptomics_pipeline.paf2blast6 import paf2blast6

def minimap2_reads(args: argparse.Namespace):
    start = time.time()
    dirpath = args.dirpath
    bigger_1 = dirpath + "/preprocessing/unassembled_reads_longer_fwd.fq"
    bigger_2 = dirpath + "/preprocessing/unassembled_reads_longer_rev.fq"
    minimap2_long_reads_out = dirpath + "/alignments/minimap2_lr_out.paf"
    tempdir = dirpath + "/alignments/minimap2_lr_alignment_outs"

    if os.path.isdir(tempdir):
        shutil.rmtree(tempdir)

    os.mkdir(tempdir)
    minimap2_command = "minimap2 -c -x sr -t " + str(args.threads) + " --split-prefix " + tempdir + "/temp" + " " + args.minimap2_index + " " + bigger_1 + " " + bigger_2 + " > " + minimap2_long_reads_out
    new_command = subprocess.run(minimap2_command, shell=True)
    check_fail("minimap2", new_command)

    end = time.time()
    print("minimap2 long read alignment against nt took: " + str(end - start))

    paf2blast6(minimap2_long_reads_out, dirpath + "/alignments")