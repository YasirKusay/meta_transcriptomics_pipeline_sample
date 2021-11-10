# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import subprocess

def merge_sams(snap_sam, diamond_sam, path):
    output_snap = open(path + "/nucl_alignments.txt", "w")
    output_diamond = open(path + "/prot_alignments.txt", "w")

    # need to firstly merge and sort the 2 files
    snap_diamond_combined_file = path + "/snap_diamond_combined_file"
    command_1 = 'egrep -v "^@" ' + snap_sam + " > " + snap_diamond_combined_file
    command = subprocess.run(command_1, shell=True) 
    command_2 = 'egrep -v "^@" ' + diamond_sam + " >> " + snap_diamond_combined_file
    command = subprocess.run(command_2, shell=True) 

    # now lets sort based on read/contig id
    snap_diamond_sorted_file = path + "/snap_diamond_sorted_file"
    sort_command = "sort -k1 " + snap_diamond_combined_file + " > " + snap_diamond_sorted_file
    command = subprocess.run(sort_command, shell=True) 

    # now we can go about selecting the best hit
    with open(snap_diamond_sorted_file, "r") as f:
        prev_query = "NULL"
        best_edit_dist = -1
        best_mapq = -1
        best_line = "NULL"
        for line in f:
            curr = line.split()
            if (curr[2] != "*"):
                assert curr[12].split(":")[0] == "NM"
                curr_line = line
                if (prev_query == curr[0]):
                    if (int(curr[12].split(":")[2]) > best_edit_dist): # checking if line has higher edit dist
                        best_line = line
                        best_edit_dist = int(curr[12].split(":")[2])  
                        best_mapq = int(curr[4])
                    elif (int(curr[12].split(":")[2]) == best_edit_dist): # checking if line has a higher mapq, if edit dists are equal
                        if (int(curr[4]) > best_mapq and int(curr[4]) != 255):
                            best_line = line
                            best_edit_dist = int(curr[12].split(":")[2])  
                            best_mapq = int(curr[4])
                    prev_query = curr[0]

                else:
                    if (best_line != "NULL"):
                        # diamond file, because it stores E-value 
                        to_print = best_line.split()
                        if to_print[15].split(":")[0] == "ZE":
                            output_diamond.write(to_print[0] + "\t" + to_print[2] + "\n")
                        else:
                            output_snap.write(to_print[0] + "\t" + to_print[2] + "\n")

                    prev_query = curr[0]
                    best_edit_dist = int(curr[12].split(":")[2])
                    best_mapq = int(curr[4])
                    best_line = line

    command = subprocess.run("rm " + snap_diamond_combined_file, shell=True)
    command = subprocess.run("rm " + snap_diamond_sorted_file, shell=True)
    
    output_snap.close()
    output_diamond.close()

    return output_snap, output_diamond