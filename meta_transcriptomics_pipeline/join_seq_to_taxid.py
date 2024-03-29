# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import os
from meta_transcriptomics_pipeline.obtain_relevant_taxids import obtain_relevant_taxids
from meta_transcriptomics_pipeline.helpers import run_shell_command

def join_accessions_taxids(hits_file, retreived_mappings, append_file):
    accessions_mapped = {}
    with open(retreived_mappings, "r") as f:         
        for line in f:
            line = line.strip()
            curr = line.split()
            accession = curr[0].split(".")[0]
            taxid = ""
            if len(curr) > 2:
                taxid = curr[2]
            else:
                # prot.accession2taxid.FULL full only has 2 columns, accession and taxid
                taxid = curr[1]
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
    run_shell_command("cat " + hits_file + " | cut -f1,2" + " | LC_COLLATE=C sort -k2  > " + hits_file_sorted)

    unique_accessions = path + "/unique_accessions.txt"
    run_shell_command("cut -f2 " + hits_file_sorted + " | uniq | egrep -v 'SNAP|Illumina|unsorted|LN:' | LC_COLLATE=C sort -k1 > " + unique_accessions)
    # for the above command, some of the nt accessions may have e.g. NR_001_name, hence cut -d'_' -f1,2 takes care of that

    retreived_mappings = path + "/retreived_mappings.txt"
    if os.path.isfile(retreived_mappings):
        os.remove(retreived_mappings)

    #for file in mapping_file:
    for curr_mapping_file in mapping_file:
        obtain_relevant_taxids(unique_accessions, curr_mapping_file, retreived_mappings)
        retreived_mappings_sorted = path + "/retreived_mappings_sorted.txt"
        if os.path.isfile(retreived_mappings) is False:
            continue
        run_shell_command("LC_COLLATE=C sort -k1 " + retreived_mappings + " > " + retreived_mappings_sorted)
        join_accessions_taxids(hits_file_sorted, retreived_mappings_sorted, append_file)
        if os.path.exists(retreived_mappings):
            os.remove(retreived_mappings)
    
def join_seq_to_taxid(nt_hits_file, nr_hits_file, nt_accession_to_taxid_mapping_files, nr_accession_to_taxid_mapping_files, contigs_reads_accessions_taxids_outfile, path):
    join(nt_hits_file, nt_accession_to_taxid_mapping_files, contigs_reads_accessions_taxids_outfile, path)
    join(nr_hits_file, nr_accession_to_taxid_mapping_files, contigs_reads_accessions_taxids_outfile, path)