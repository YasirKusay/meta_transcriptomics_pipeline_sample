import argparse
import subprocess
import time
from meta_transcriptomics_pipeline.helpers import check_fail

def run_diamond(index, in_path, out_path, threads):
    diamond_command = "diamond" + " blastx --db " + index +\
                        " --query " + in_path + " --mid-sensitive --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" +\
                        " --masking 0 -c 1 -b 6" +\
                        " --threads " + str(threads) +\
                        " --out " + out_path
    new_command = subprocess.run(diamond_command, shell=True)
    if check_fail("diamond", new_command) is True: return False

def diamond(args: argparse.Namespace):
    dirpath = args.dirpath

    if dirpath[-1] == "/":
        dirpath = dirpath[0:-1]

    nr_combined_file = dirpath + "/alignments/nr_alignments_file.tsv"
    combined_file = dirpath + "/preprocessing/combined_reads_contigs_file.fq"

    start = time.time()
    if run_diamond(args.diamond_index, combined_file, nr_combined_file, args.threads) == False: return None
    end = time.time()
    print("contig alignment against nr took: " + str(end - start))