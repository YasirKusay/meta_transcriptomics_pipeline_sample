import argparse
import subprocess
import time
import os
from helpers import check_fail
from paf2blast6 import paf2blast6

def minimap2_contigs(args: argparse.Namespace):
    start = time.time()
    dirpath = args.dirpath
    minimap2_path = "minimap2"
    contig_path = dirpath + "/megahit_out"
    minimap2_contig_out = dirpath + "/minimap2_contig_out.paf"
    new_contigs = contig_path + "/final_contigs.fq"
    tempdir = dirpath + "/tempdir1"
    os.mkdir(tempdir)
    minimap2_command = minimap2_path + " -c -t " + str(args.threads) + " --split-prefix " + tempdir + "/temp" + " " + args.minimap2_index + " " + new_contigs + " > " + minimap2_contig_out
    new_command = subprocess.run(minimap2_command, shell=True)
    if check_fail(minimap2_path, new_command, []) is True: return None
    end = time.time()
    print("minimap2 alignment 1 against nt took: " + str(end - start))
    
    paf2blast6(minimap2_contig_out, dirpath)
    paf2blast_out = dirpath + "/minimap2_contig_out_frompaf.m8"