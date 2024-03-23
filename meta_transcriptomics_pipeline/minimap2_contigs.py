import argparse
import logging
import time
import os
import shutil
from meta_transcriptomics_pipeline.helpers import run_shell_command, check_megahit_success
from meta_transcriptomics_pipeline.paf2blast6 import paf2blast6

log = logging.getLogger(__name__)

def minimap2_contigs(args: argparse.Namespace):
    start = time.time()
    dirpath = args.dirpath
    contig_path = dirpath + "/megahit_out"
    minimap2_contig_out = dirpath + "/alignments/minimap2_contig_out.paf"
    new_contigs = contig_path + "/final_contigs.fq"
    tempdir = dirpath + "/alignments/minimap2_contigs_alignment_outs"

    if os.path.isdir(tempdir):
        shutil.rmtree(tempdir)

    # I want to stop this step without stopping the pipeline itself if megahit failed
    # I did it this way as a temporary solution as I would otherwise need to rewrite
    # The submission scripts (it is too difficult to to update given there are multiple)
    # Copies of it
    if check_megahit_success(contig_path + "/log") == False:
        log.error("Skipping minimap2_contigs as megahit had failed.")
        return

    os.mkdir(tempdir)

    minimap2_command = "minimap2 -c -t " + str(args.threads) + " --split-prefix " + tempdir + "/temp" + " " + args.minimap2_index + " " + new_contigs + " > " + minimap2_contig_out
    run_shell_command(minimap2_command)
    end = time.time()
    log.info("minimap2 contig alignment against nt took: " + str(end - start))
    
    paf2blast6(minimap2_contig_out, dirpath + "/alignments")