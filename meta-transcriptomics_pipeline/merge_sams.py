# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import sys
import os
import subprocess
from helpers import check_fail

def merge_fasta(snap_sam, diamond_sam, path, seq_type):
    # lets first merge the 2 files together
    out = path + "/merged_temp" + seq_type + ".sam" 
    merge_command = "samtools merge " + out + " " + snap_sam + " " + diamond_sam 
    new_command = subprocess.run(merge_command, shell=True)
    if check_fail("samtools merge", new_command, []) is True: return None

    final_file = path + "/merged_" + seq_type + ".sam" 
    output = open(final_file, "w")

    # now we can go about selecting the best hit

    with open(out, "r") as f:
        prev_query = "NULL"
        best_query = "NULL"
        best_edit_dist = -1
        best_mapq = -1
        best_line = "NULL"
        for line in f:
            curr = line.split()
            assert curr[12].split(":")[0] == "NM"
            if (curr[0][0]!="@"):
                curr_line = line
                if (prev_query == curr[0]):
                    if (int(curr[12].split(":")[2]) > best_edit_dist): # checking if line has higher edit dist
                        best_query = curr[0]
                        best_line = line
                        best_edit_dist = int(curr[12].split(":")[2])  
                    elif (int(curr[12].split(":")[2]) == best_edit_dist): # checking if line has a higher mapq, if edit dists are equal
                        if (int(curr[4]) > best_mapq and int(curr[4]) != 255):
                            best_query = curr[0]
                            best_line = line
                            best_edit_dist = int(curr[12].split(":")[2])  
                            best_mapq = int(curr[4])
                    
                    pre_query = curr[0]

                else:
                    if (best_line != "NULL"):
                        output.write(line)
                        
                    prev_query = curr[0]
                    best_query = prev_query
                    best_edit_dist = int(curr[12].split(":")[2])
                    best_mapq = int(curr[4])
                    best_line = line
                    

                # we are in the alignment section

            elif (curr[0][0]=="@"):
                output.write(line)
                # we are in the header section


    output.close()

    # returning merged file name
    return out 