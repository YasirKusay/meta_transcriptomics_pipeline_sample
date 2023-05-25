import argparse
import time
from meta_transcriptomics_pipeline.helpers import run_shell_command

# explanation of how diamond flags work
# if you are experiencing memory issues, its better to raise the chunks (-c) rather than lower the block size (-b)
# it should not impact the performance too much according to this:
# https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.3176/MediaObjects/41592_2015_BFnmeth3176_MOESM5_ESM.pdf
# according to the diamond paper -B refers to the max amount of the index/query sequence letters to load and compare at a time
# -c refers to how many chunks we want to divide an individual seed when processing it 
# also the paper (towards the end) gives a formula for estimating the memory usage
# 2(B + 8 x B/C + const). Hence, for -b 6 and -c 1, the estimated memory usage is 2(6e9 + 8*6e9) = 108 GB expected mem usage

def run_diamond(index, in_path, out_path, threads):
    diamond_command = "diamond" + " blastx --db " + index +\
                        " --query " + in_path + " --mid-sensitive --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" +\
                        " --masking 0 -c 1 -b 6" +\
                        " --threads " + str(threads) +\
                        " --out " + out_path
    
    run_shell_command(diamond_command)

def diamond(args: argparse.Namespace):
    dirpath = args.dirpath

    if dirpath[-1] == "/":
        dirpath = dirpath[0:-1]

    nr_combined_file = dirpath + "/alignments/nr_alignments_file.tsv"
    combined_file = dirpath + "/preprocessing/combined_reads_contigs_file.fq"

    start = time.time()
    run_diamond(args.diamond_index, combined_file, nr_combined_file, args.threads)
    end = time.time()
    print("contig alignment against nr took: " + str(end - start))