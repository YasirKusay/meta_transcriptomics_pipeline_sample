# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import subprocess
import os

def merge_sams(snap_sam, diamond_sam, path, snap_out = None, diamond_out = None):

    if snap_out is None:
        snap_out = path + "/nucl_alignments_reads.txt"
    if diamond_out is None:
        diamond_out = path + "/prot_alignments_reads.txt"

    if os.path.isfile(snap_out):
        os.remove(snap_out)
    if os.path.isfile(diamond_out):
        os.remove(diamond_out) 

    output_snap = open(snap_out, "w")
    output_diamond = open(diamond_out, "w")

    # need to firstly merge and sort the 2 files
    snap_diamond_combined_file = path + "/snap_diamond_combined_file"
    command_1 = 'egrep -v "^@" ' + snap_sam + " > " + snap_diamond_combined_file
    command = subprocess.run(command_1, shell=True) 
    command_2 = 'egrep -v "^@" ' + diamond_sam + " >> " + snap_diamond_combined_file
    command = subprocess.run(command_2, shell=True) 

    # now lets sort based on read/contig id
    snap_diamond_sorted_file = path + "/snap_diamond_sorted_file"
    sort_command = "LC_COLLATE=C sort -k1 " + snap_diamond_combined_file + " > " + snap_diamond_sorted_file
    command = subprocess.run(sort_command, shell=True) 

    # now we can go about selecting the best hit
    with open(snap_diamond_sorted_file, "r") as f:
        prev_query = "NULL"
        best_edit_dist = -1
        best_mapq = -1
        best_line = "NULL"
        for line in f:
            curr = line.split("\t")
            print(curr)
            if (curr[2] != "*"):
                edit_dist_inc = 11
                while (edit_dist_inc < len(curr)):
                    if curr[edit_dist_inc].split(":")[0] == "NM":
                        break
                    edit_dist_inc += 1
                curr_line = line
                curr_edit_dist = -1
                if edit_dist_inc < len(curr):
                    curr_edit_dist = int(curr[edit_dist_inc].split(":")[2])
                if (prev_query == curr[0]):
                    if (curr_edit_dist > best_edit_dist): # checking if line has higher edit dist
                        best_line = line
                        best_edit_dist = curr_edit_dist
                        best_mapq = int(curr[4])
                    elif (curr_edit_dist == best_edit_dist): # checking if line has a higher mapq, if edit dists are equal
                        if (best_mapq != 255):
                            if (int(curr[4]) > best_mapq):
                                best_line = line
                                best_edit_dist = curr_edit_dist
                                best_mapq = int(curr[4])
                        else: 
                            if (int(curr[4]) != 255):
                                best_line = line
                                best_edit_dist = curr_edit_dist
                                best_mapq = int(curr[4])
                    prev_query = curr[0]

                else:
                    if (best_line != "NULL"):
                        to_print = best_line.split("\t")
                        accession = to_print[2]
                        full_accession = accession.split("_")
                        actual_accession = full_accession[:2]
                        print_accession = "_".join(actual_accession)

                        # diamond file, because it stores E-value 
                        if (15 < len(curr)):
                            if to_print[15].split(":")[0] == "ZE":
                                output_diamond.write(to_print[0] + "\t" + print_accession + "\n")
                            else:
                                output_snap.write(to_print[0] + "\t" + print_accession + "\n")
                        else: 
                            output_snap.write(to_print[0] + "\t" + print_accession + "\n")

                    prev_query = curr[0]
                    best_edit_dist = curr_edit_dist
                    best_mapq = int(curr[4])
                    best_line = line
        
        if (best_line != "NULL"):
            to_print = best_line.split("\t")
            accession = to_print[2]
            full_accession = accession.split("_")
            actual_accession = full_accession[:2]
            print_accession = "_".join(actual_accession)

            # diamond file, because it stores E-value 
            if (15 < len(curr)):
                if to_print[15].split(":")[0] == "ZE":
                    output_diamond.write(to_print[0] + "\t" + print_accession + "\n")
                else:
                    output_snap.write(to_print[0] + "\t" + print_accession + "\n")
            else: 
                output_snap.write(to_print[0] + "\t" + print_accession + "\n")

    #command = subprocess.run("rm " + snap_diamond_combined_file, shell=True)
    #command = subprocess.run("rm " + snap_diamond_sorted_file, shell=True)
    
    output_snap.close()
    output_diamond.close()

    return snap_out, diamond_out
