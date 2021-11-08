# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import subprocess
from helpers import check_fail

def obtain_relevant_taxids(accession_file, mapping_file, write_file):

    mf = open(accession_file, "r")
    curr_accession = mf.readline()
    wf = open(write_file, "w")
    
    with open(mapping_file, "r") as cf:
        for line in cf:
            curr = line.split()
            if (curr[0] == curr_accession):
                wf.write(line)
                curr_accession = mf.readline()

    mf.close()
    wf.close()
            

def join(contig_file, read_file, mapping_file, join_file, path):

    # lets firstly join the contig and read files 
    combined_file_temp = path + "/combined_file_temp.txt"
    combined_file_sorted = path + "/combined_file_sorted.txt"
    subprocess.run("cat " + contig_file + " > " + combined_file_temp, shell=True)
    subprocess.run("cat " + read_file + " >> " + combined_file_temp, shell=True)
    subprocess.run("sort -k2 " + combined_file_temp + " > " + combined_file_sorted, shell=True)

    unique_accessions = path + "/unique_accessions.txt"
    subprocess.run("cut -f2 | uniq " + combined_file_sorted + " > " + unique_accessions, shell=True)

    relevant_taxids = path + "/relevant_taxids.txt"
    obtain_relevant_taxids(unique_accessions, mapping_file, relevant_taxids)

    command = "join -1 2 -2 1 -o \"1.1, 1.2, 2.1\" <(sort -k2 " + combined_file_sorted + ") <(sort -k1 " + relevant_taxids + ") >> " + join_file
    new_command = subprocess.run(command, shell=True)
    if check_fail("join", new_command, []) is True: return None
    
def join_taxid_contigs(output_contig_snap, output_read_snap, output_contig_diamond, output_read_diamond, nucl_map, prot_map, final_out, path):
    if join(output_contig_snap, output_read_snap, nucl_map, final_out, path) == None or join(output_contig_diamond, output_read_diamond, prot_map, final_out, path) == None: 
        return False

    return True
    