# file takes the 2 blast files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by e value (1) and bit scores (2)

import subprocess
import os

def merge_contigs(snap_sam, diamond_sam, path):
    if os.path.isfile(path + "/nucl_alignments_contigs.txt"):
        os.remove(path + "/nucl_alignments_contigs.txt")
    if os.path.isfile(path + "/prot_alignments_contigs.txt"):
        os.remove(path + "/prot_alignments_contigs.txt")    

    output_snap = open(path + "/nucl_alignments_contigs.txt", "w")
    output_diamond = open(path + "/prot_alignments_contigs.txt", "w")

    # need to firstly merge and sort the 2 files
    temp_sam = path + "/temp_sam"
    add_n = "awk '{print $0, \"\tN\"}' " + snap_sam + " > " + temp_sam
    command = subprocess.run(add_n, shell=True)
    temp_diamond = path + "/temp_diamond"
    add_p = "awk '{print $0, \"\tP\"}' " + diamond_sam + " > " + temp_diamond
    command = subprocess.run(add_p, shell=True)

    snap_diamond_combined_file = path + "/snap_diamond_combined_file"
    command_1 = 'cat ' + temp_sam + " > " + snap_diamond_combined_file
    command = subprocess.run(command_1, shell=True) 
    command_2 = 'cat ' + temp_diamond + " >> " + snap_diamond_combined_file
    command = subprocess.run(command_2, shell=True)

    # N means Nucleotide, P means protein

    # now lets sort based on read/contig id
    snap_diamond_sorted_file = path + "/snap_diamond_sorted_file"
    sort_command = "LC_COLLATE=C sort -k1 " + snap_diamond_combined_file + " > " + snap_diamond_sorted_file
    command = subprocess.run(sort_command, shell=True) 

    # now we can go about selecting the best hit
    with open(snap_diamond_sorted_file, "r") as f:
        prev_query = "NULL"
        best_e_value = -1
        best_bitscore = -1
        best_line = "NULL"
        for line in f:
            curr = line.split("\t")
            curr_line = line
            if (prev_query == curr[0]):
                if (float(curr[10]) > best_e_value): # e value
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
                    accession = to_print[1]
                    full_accession = accession.split("_")
                    actual_accession = full_accession[:2]
                    print_accession = "_".join(actual_accession)
                    # diamond file, because it stores E-value
                    print(to_print[12]) 
                    if to_print[12].strip() == "P":
                        output_diamond.write(to_print[0] + "\t" + print_accession + "\t" + to_print[10] + "\t" + to_print[11] + "\t" + to_print[2] + "\n")
                    else:
                        output_snap.write(to_print[0] + "\t" + print_accession + "\t" + to_print[10] + "\t" + to_print[11] + "\t" + to_print[2] + "\n")

                prev_query = curr[0]
                best_e_value = float(curr[10])
                best_bitscore = float(curr[11])
                best_line = line

        if (best_line != "NULL"):
            to_print = best_line.split("\t")
            accession = to_print[1]
            full_accession = accession.split("_")
            actual_accession = full_accession[:2]
            print_accession = "_".join(actual_accession)
            # diamond file, because it stores E-value
            print(to_print[12]) 
            if to_print[12].strip() == "P":
                output_diamond.write(to_print[0] + "\t" + print_accession + "\t" + to_print[10] + "\t" + to_print[11] + "\t" + to_print[2] + "\n")
            else:
                output_snap.write(to_print[0] + "\t" + print_accession + "\t" + to_print[10] + "\t" + to_print[11] + "\t" + to_print[2] + "\n")
                
    #command = subprocess.run("rm " + snap_diamond_combined_file, shell=True)
    #command = subprocess.run("rm " + snap_diamond_sorted_file, shell=True)
    
    output_snap.close()
    output_diamond.close()

    return path + "/nucl_alignments_contigs.txt", path + "/prot_alignments_contigs.txt"
