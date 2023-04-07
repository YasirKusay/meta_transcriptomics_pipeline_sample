import argparse
import subprocess
import time
import os
from meta_transcriptomics_pipeline.helpers import check_fail

def blast_sr(args: argparse.Namespace):
    dirpath = args.dirpath
    if dirpath[-1] == "/":
        dirpath = dirpath[0:-1]

    blast_path = "blastn -task megablast"

    nt_alignments_file = dirpath + "/alignments/nt_alignments_file.tsv"
    if os.path.isfile(nt_alignments_file):
        os.remove(nt_alignments_file)

    combined_file_fa = dirpath + "/preprocessing/combined_reads_contigs_file.fa"
    
    start = time.time()
    blast_command = blast_path + " -query " + combined_file_fa + " -db " + args.blast_ncbi_nt_database + " -out " + nt_alignments_file + " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen\" -max_target_seqs 10 -num_threads " + str(args.threads)
    new_command = subprocess.run(blast_command, shell=True)
    if check_fail(blast_path, new_command) is True: return False  

    end = time.time()
    print("blast alignment against nt took: " + str(end - start))