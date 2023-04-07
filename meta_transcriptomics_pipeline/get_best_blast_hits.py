# file takes the 2 blast files, one from blast (and minimap) and the other for diamond and merges them
# it finds the best hits possible
# for each hit, it will place it in a file depending on whether the hit originated from the nt
# mapping file or the nr mapping file
# best hits are decided by e value (1) and bit scores (2)

import subprocess
import os

def get_best_blast_hits(nt_alignments_file, nr_alignments_file, path, best_nt_outfile, best_nr_outfile):

    best_nt_output = open(best_nt_outfile, "w")
    best_nr_output = open(best_nr_outfile, "w")

    # need to firstly merge and sort the 2 files
    temp_nt = path + "/temp_nt"
    add_nt = "awk '{print $0, \"\tN\"}' " + nt_alignments_file + " > " + temp_nt
    command = subprocess.run(add_nt, shell=True)

    temp_nr = path + "/temp_nr"
    add_nr = "awk '{print $0, \"\tP\"}' " + nr_alignments_file + " > " + temp_nr
    command = subprocess.run(add_nr, shell=True)

    nt_nr_alignments_combined_file = path + "/nt_nr_alignments_combined"
    command_1 = 'cat ' + temp_nt + " > " + nt_nr_alignments_combined_file
    command = subprocess.run(command_1, shell=True) 
    command_2 = 'cat ' + temp_nr + " >> " + nt_nr_alignments_combined_file
    command = subprocess.run(command_2, shell=True)

    os.remove(temp_nt)
    os.remove(temp_nr)

    # N means Nucleotide, P means protein

    # now lets sort based on read/contig id
    nt_nr_alignments_combined_file_sorted = path + "/nt_nr_alignments_combined_sorted"
    sort_command = "LC_COLLATE=C sort -k1 " + nt_nr_alignments_combined_file + " > " + nt_nr_alignments_combined_file_sorted
    command = subprocess.run(sort_command, shell=True) 

    # now we can go about selecting the best hit
    with open(nt_nr_alignments_combined_file_sorted, "r") as f:
        prev_query = "NULL"
        best_e_value = -1
        best_bitscore = -1
        best_line = "NULL"
        for line in f:
            curr = line.split("\t")
            if (prev_query == curr[0]):
                if (float(curr[10]) < best_e_value): # e value
                    best_line = line
                    best_e_value = float(curr[10])
                    best_bitscore = float(curr[11])
                elif (float(curr[10]) == best_e_value): # checking if line has a higher e value, if bit scores
                    if (float(curr[11]) > best_bitscore ):
                        best_line = line
                        best_e_value = float(curr[10])
                        best_bitscore = float(curr[11])
                prev_query = curr[0]

            else:
                if (best_line != "NULL"):
                    to_print = best_line.split("\t")
                    read = to_print[0]
                    accession = to_print[1].split(".")[0] # we only care about the accession, not the version
                    percent_id = to_print[2]
                    e_value = to_print[10]
                    bitscore = to_print[11].strip()
                    qlen = to_print[12].strip()

                    if to_print[13].strip() == "P": # diamond file
                        best_nr_output.write(read + "\t" + accession + "\t" + e_value + "\t" + bitscore + "\t" + percent_id + "\t" + qlen + "\n")
                    else:
                        best_nt_output.write(read + "\t" + accession + "\t" + e_value + "\t" + bitscore + "\t" + percent_id + "\t" + qlen + "\n")

                prev_query = curr[0]
                best_e_value = float(curr[10])
                best_bitscore = float(curr[11])
                best_line = line

        if (best_line != "NULL"):
            to_print = best_line.split("\t")
            read = to_print[0]
            accession = to_print[1].split(".")[0]
            percent_id = to_print[2]
            e_value = to_print[10]
            bitscore = to_print[11].strip()
            qlen = to_print[12].strip()

            if to_print[13].strip() == "P": # diamond file
                best_nr_output.write(read + "\t" + accession + "\t" + e_value + "\t" + bitscore + "\t" + percent_id + "\t" + qlen + "\n")
            else:
                best_nt_output.write(read + "\t" + accession + "\t" + e_value + "\t" + bitscore + "\t" + percent_id + "\t" + qlen + "\n")
    
    best_nt_output.close()
    best_nr_output.close()
