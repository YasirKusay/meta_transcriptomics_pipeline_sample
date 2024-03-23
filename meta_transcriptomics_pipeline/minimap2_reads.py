import argparse
import logging
import time
import os
import shutil
from meta_transcriptomics_pipeline.helpers import run_shell_command
from meta_transcriptomics_pipeline.paf2blast6 import paf2blast6

log = logging.getLogger(__name__)

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
    run_shell_command(minimap2_command)

    end = time.time()
    log.info("minimap2 long read alignment against nt took: " + str(end - start))

    paf2blast6(minimap2_long_reads_out, dirpath + "/alignments")