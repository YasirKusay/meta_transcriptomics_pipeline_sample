import argparse
import subprocess
import time
from meta_transcriptomics_pipeline.helpers import check_fail

def blast_sr(args: argparse.Namespace):
    dirpath = args.dirpath

    smaller_1 = dirpath + "/smaller_1.fq"
    smaller_2 = dirpath + "/smaller_2.fq"
    # now lets blast the small reads, (merge them first)

    seqtk_path = "seqtk"
    
    merged_short_pe = dirpath + "/merged_short_reads.fq"
    merge_command = "seqtk mergepe " + smaller_1 + " " + smaller_2 + " > " + merged_short_pe
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("seqtk mergepe", new_command, []) is True: return False 

    merged_short_pe_fa = dirpath + "/merged_short_reads.fa"
    seqtk_command = seqtk_path + " seq -a " + merged_short_pe + " > " + merged_short_pe_fa
    new_command = subprocess.run(seqtk_command, shell=True)
    if check_fail("seqtk", new_command, []) is True: return False 

    start = time.time()

    blast_path = "blastn -task megablast"
    blast_short_out = dirpath + "/blast_short_out"
    blast_command = blast_path + " -query " + merged_short_pe_fa + " -db nt " + " -out " + blast_short_out + " -outfmt 6 -max_target_seqs 20 -num_threads " + str(args.threads)
    new_command = subprocess.run(blast_command, shell=True)
    if check_fail(blast_path, new_command, []) is True: return False  

    end = time.time()
    print("blast alignment against nt took: " + str(end - start))

    nt_combined_file = dirpath + "/nt_combined_file"

    new_command = subprocess.run("cat " + blast_short_out + " >> " + nt_combined_file, shell=True)