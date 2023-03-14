# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import subprocess
import os
from meta_transcriptomics_pipeline.obtain_relevant_taxids import obtain_relevant_taxids

def join_accessions_taxids(hits_file, retreived_mappings, append_file):
    accessions_mapped = {}
    with open(retreived_mappings, "r") as f:         
        for line in f:
            curr = line.split()
            accession = curr[0]
            taxid = curr[2]
            accessions_mapped[accession] = taxid

    wf = open(append_file, "a")
    with open(hits_file, "r") as f:
        for line in f:
            curr = line.split()            
            read = curr[0]
            full_accession = curr[1]
            accession = full_accession.split(".")[0]
            if (accession in accessions_mapped.keys()):
                wf.write(read + "\t" + accession + "\t" + accessions_mapped[accession] + "\n")

    wf.close()

def join(hits_file, mapping_file, append_file, path):

    # lets firstly join the contig and read files 
    hits_file_sorted = path + "/hits_file_sorted.txt"
    new_command = subprocess.run("cat " + hits_file + " | cut -f1,2" + " | LC_COLLATE=C sort -k2  > " + hits_file_sorted, shell=True)

    unique_accessions = path + "/unique_accessions.txt"
    new_command = subprocess.run("cut -f2 " + hits_file_sorted + " | uniq | egrep -v 'SNAP|Illumina|unsorted|LN:' | LC_COLLATE=C sort -k1 > " + unique_accessions, shell=True)
    # for the above command, some of the nt accessions may have e.g. NR_001_name, hence cut -d'_' -f1,2 takes care of that

    retreived_mappings = path + "/retreived_mappings.txt"
    if os.path.isfile(retreived_mappings):
        os.remove(retreived_mappings)

    #for file in mapping_file:
    for curr_mapping_file in mapping_file:
        obtain_relevant_taxids(unique_accessions, curr_mapping_file, retreived_mappings)
        retreived_mappings_sorted = path + "/retreived_mappings_sorted.txt"
        new_command = subprocess.run("LC_COLLATE=C sort -k1 " + retreived_mappings + " > " + retreived_mappings_sorted, shell = True)
        join_accessions_taxids(hits_file_sorted, retreived_mappings_sorted, append_file)
        if os.path.exists(retreived_mappings):
            os.remove(retreived_mappings)

    return True
    
def join_taxid_contigs(nt_hits_file, nr_hits_file, nt_accession_to_taxid_mapping_files, nr_accession_to_taxid_mapping_files, contigs_reads_accessions_taxids_outfile, path):
    if join(nt_hits_file, nt_accession_to_taxid_mapping_files, contigs_reads_accessions_taxids_outfile, path) is False or join(nr_hits_file, nr_accession_to_taxid_mapping_files, contigs_reads_accessions_taxids_outfile, path) is False: 
        return False

    return True