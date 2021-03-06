# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import subprocess
import os
from helpers import check_fail

def obtain_relevant_taxids(accession_file, mapping_file, write_file):

    mf = open(accession_file, "r")
    curr_accession = mf.readline().strip()
    wf = open(write_file, "a")

    while True:
        if (curr_accession is None or curr_accession == "" or curr_accession == "*"):
            curr_accession = mf.readline().rstrip()
        else:
            break
    
    with open(mapping_file, "r") as cf:
        for line in cf:
            curr = line.split()
            if (curr[0] == curr_accession):
                wf.write(line)
                curr_accession = mf.readline().strip()
                if (curr_accession is None or curr_accession == "" or curr_accession == "\n" or curr_accession == "\r\n" or len(curr_accession) == 0):
                    break

            if (curr[0][0] > curr_accession[0]): # check if first character is greater in taxlist
                curr_accession = mf.readline().strip()
                if (curr_accession is None or curr_accession == "" or curr_accession == "\n" or curr_accession == "\r\n" or len(curr_accession) == 0):
                    break

    mf.close()
    wf.close()
            

def join(contig_file, read_file, mapping_file, join_file, path):

    # lets firstly join the contig and read files 
    combined_file_temp = path + "/combined_file_temp.txt"
    combined_file_sorted = path + "/combined_file_sorted.txt"
    new_command = subprocess.run("cat " + contig_file + " | cut -f1,2" + " > " + combined_file_temp, shell=True)
    new_command = subprocess.run("cat " + read_file + " | cut -f1,2" + " >> " + combined_file_temp, shell=True)
    new_command = subprocess.run("sort -k2 " + combined_file_temp + " > " + combined_file_sorted, shell=True)

    unique_accessions = path + "/unique_accessions.txt"
    new_command = subprocess.run("cut -f2 " + combined_file_sorted + " | uniq | egrep -v 'SNAP|Illumina|unsorted|LN:' | cut -d'_' -f1,2  > " + unique_accessions, shell=True)
    # for the above command, some of the nt accessions may have e.g. NR_001_name, hence cut -d'_' -f1,2 takes care of that

    relevant_taxids = path + "/relevant_taxids.txt"
    obtain_relevant_taxids(unique_accessions, mapping_file, relevant_taxids)

    combined_file_sorted_sorted = path + "/combined_file_sorted_sorted"
    new_command = subprocess.run("sort -k2 " + combined_file_sorted + " > " + combined_file_sorted_sorted, shell = True)
    relevant_taxids_sorted = path + "/relevant_taxids_sorted"
    new_command = subprocess.run("sort -k1 " + relevant_taxids + " > " + relevant_taxids_sorted, shell = True)

    command = ""
    if (os.path.isfile(join_file)):
        command = "join -1 2 -2 1 -o \"1.1 1.2 2.2\" " + combined_file_sorted_sorted + " " + relevant_taxids_sorted +  " >> " + join_file
    else:
        command = "join -1 2 -2 1 -o \"1.1 1.2 2.2\" " + combined_file_sorted_sorted + " " + relevant_taxids_sorted +  " > " + join_file
    new_command = subprocess.run(command, shell=True)
    #if check_fail("join", new_command, []) is True: return None
    return True
    
def join_taxid_contigs(output_contig_snap, output_read_snap, output_contig_diamond, output_read_diamond, nucl_map, prot_map, final_out, path):
    if join(output_contig_snap, output_read_snap, nucl_map, final_out, path) is False or join(output_contig_diamond, output_read_diamond, prot_map, final_out, path) is False: 
        return False

    return True
    