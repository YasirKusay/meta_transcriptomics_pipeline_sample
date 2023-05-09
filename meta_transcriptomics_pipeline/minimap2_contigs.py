import argparse
import subprocess
import time
import os
from meta_transcriptomics_pipeline.helpers import check_fail
from meta_transcriptomics_pipeline.paf2blast6 import paf2blast6

def minimap2_contigs(args: argparse.Namespace):
    start = time.time()
    dirpath = args.dirpath
    contig_path = dirpath + "/megahit_out"
    minimap2_contig_out = dirpath + "/alignments/minimap2_contig_out.paf"
    new_contigs = contig_path + "/final_contigs.fq"
    tempdir = dirpath + "/alignments/minimap2_contigs_alignment_outs"
    os.mkdir(tempdir)
    minimap2_command = "minimap2 -c -t " + str(args.threads) + " --split-prefix " + tempdir + "/temp" + " " + args.minimap2_index + " " + new_contigs + " > " + minimap2_contig_out
    new_command = subprocess.run(minimap2_command, shell=True)
    if check_fail("minimap2", new_command) is True: return None
    end = time.time()
    print("minimap2 contig alignment against nt took: " + str(end - start))
    
    paf2blast6(minimap2_contig_out, dirpath + "/alignments")