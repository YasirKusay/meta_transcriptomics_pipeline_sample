# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import subprocess
import os
from meta_transcriptomics_pipeline.helpers import check_fail
from meta_transcriptomics_pipeline.obtain_relevant_taxids import obtain_relevant_taxids

def join_accessions_taxids(combined_file_sorted, relevant_taxids, out):
    accessions_mapped = {}
    with open(relevant_taxids, "r") as f:         
        for line in f:
            curr = line.split()
            accession = curr[0]
            taxid = curr[1]
            accessions_mapped[accession] = taxid

    wf = open(out, "a")
    with open(combined_file_sorted, "r") as f:
        for line in f:
            curr = line.split()            
            read = curr[0]
            accession = curr[1]
            if (accession in accessions_mapped.keys()):
                wf.write(read + "\t" + accession + "\t" + accessions_mapped[accession] + "\n")

    wf.close()

def join(contig_file, read_file, mapping_file, join_file, path):

    # lets firstly join the contig and read files 
    combined_file_temp = path + "/combined_file_temp.txt"
    combined_file_sorted = path + "/combined_file_sorted.txt"
    new_command = subprocess.run("cat " + contig_file + " | cut -f1,2" + " > " + combined_file_temp, shell=True)
    new_command = subprocess.run("cat " + read_file + " | cut -f1,2" + " >> " + combined_file_temp, shell=True)
    new_command = subprocess.run("LC_COLLATE=C sort -k2 " + combined_file_temp + " > " + combined_file_sorted, shell=True)

    unique_accessions = path + "/unique_accessions.txt"
    new_command = subprocess.run("cut -f2 " + combined_file_sorted + " | uniq | egrep -v 'SNAP|Illumina|unsorted|LN:' | LC_COLLATE=C sort -k1  > " + unique_accessions, shell=True)
    # for the above command, some of the nt accessions may have e.g. NR_001_name, hence cut -d'_' -f1,2 takes care of that
    
    relevant_taxids = path + "/relevant_taxids.txt"
    obtain_relevant_taxids(unique_accessions, mapping_file, relevant_taxids)

    combined_file_sorted_sorted = path + "/combined_file_sorted_sorted"
    new_command = subprocess.run("LC_COLLATE=C sort -k2 " + combined_file_sorted + " > " + combined_file_sorted_sorted, shell = True)
    relevant_taxids_sorted = path + "/relevant_taxids_sorted"
    new_command = subprocess.run("LC_COLLATE=C sort -k1 " + relevant_taxids + " > " + relevant_taxids_sorted, shell = True)

    join_accessions_taxids(combined_file_sorted_sorted, relevant_taxids_sorted, join_file)
    return True
    
def join_taxid_contigs(output_contig_snap, output_read_snap, output_contig_diamond, output_read_diamond, nucl_map, prot_map, final_out, path):
    relevant_taxids = path + "/relevant_taxids.txt"
    if os.path.isfile(relevant_taxids):
        os.remove(relevant_taxids)
    if join(output_contig_snap, output_read_snap, nucl_map, final_out, path) is False or join(output_contig_diamond, output_read_diamond, prot_map, final_out, path) is False: 
        return False

    return True
    
