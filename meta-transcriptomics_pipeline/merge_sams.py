# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import sys
import os
import subprocess

def merge_fasta(snap_sam, diamond_sam, path):
    output_snap = open(path + "/nucl_alignments.txt", "w")
    output_diamond = open(path + "/prot_alignments.txt", "w")
    best_hits = {}

    # now we can go about selecting the best hit
    for out in [snap_sam, diamond_sam]:
        with open(out, "r") as f:
            for line in f:
                curr = line.split()
                assert curr[12].split(":")[0] == "NM"
                if (curr[0][0]!="@"):
                    curr_line = line
                    if (curr[0] in best_hits):
                        prev_best = best_hits[curr[0]][1]
                        prev = prev_best.split()
                        if (int(curr[12].split(":")[2]) > int(prev[12].split(":")[2])): # checking if line has higher edit dist
                            best_hits[curr[0]] = [out, line]
                        elif (int(curr[12].split(":")[2]) == int(prev[12].split(":")[2])): # checking if line has a higher mapq, if edit dists are equal
                            if (int(curr[4]) > int(prev[4]) and int(curr[4]) != 255):
                                best_hits[curr[0]] = [out, line]
                        
                        pre_query = curr[0]

                    else:   
                        best_hits[curr[0]] = [out, line]

    
    for key in best_hits:
        line = best_hits[key][1]
        curr = line.split()
        if best_hits[key][0] == snap_sam:
            output_snap.write(key + '\t' + curr[2])
        else: 
            output_diamond.write(key + '\t' + curr[2])
    
    output_snap.close()
    output_diamond.close()

    return output_snap, output_diamond