import argparse
import subprocess
import time
import os
from meta_transcriptomics_pipeline.helpers import check_fail

def blast_sr(args: argparse.Namespace):
    dirpath = args.dirpath

    start = time.time()

    blast_path = "blastn -task megablast"
    blast_short_out = dirpath + "/blast_short_out"
    combined_file_fa = dirpath + "/combined_file.fa"
    
    blast_command = blast_path + " -query " + combined_file_fa + " -db " + args.blast_ncbi_nt_database + " -out " + blast_short_out + " -outfmt 6 -max_target_seqs 20 -num_threads " + str(args.threads)
    new_command = subprocess.run(blast_command, shell=True)
    if check_fail(blast_path, new_command, []) is True: return False  

    end = time.time()
    print("blast alignment against nt took: " + str(end - start))

    nt_combined_file = dirpath + "/nt_combined_file"
    if os.path.isfile(nt_combined_file):
        os.remove(nt_combined_file)

    new_command = subprocess.run("cat " + blast_short_out + " >> " + nt_combined_file, shell=True)