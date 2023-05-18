import os
from meta_transcriptomics_pipeline.helpers import run_shell_command

# takes a sam file
# returns a file matching a read with its best matching hit (if a match exists)
# best hits are decided by lowest edit distance (1) and highest mapq (2)
# the path places intermediate outputs
def get_best_sam_hits(snap_file, best_read_hit_file, path):

    if os.path.isfile(best_read_hit_file):
        os.remove(best_read_hit_file)

    best_hits_output = open(best_read_hit_file, "w")

    # need to firstly merge and sort the 2 files based on read id
    snap_sorted_file = path + "/snap_sorted_file"
    combine_and_sort_snap_command = 'egrep -v "^@" ' + snap_file + " | LC_COLLATE=C sort -k1 > " + snap_sorted_file
    run_shell_command(combine_and_sort_snap_command) 

    # now we can go about selecting the best hit
    # column 1 is the read name
    # column 3 is the item the read mapped onto (* if no mapping exists)
    # column 5 is the mapping quality
    # edit_distance is the column beginning with NM: (somewhere after column 11)
    with open(snap_sorted_file, "r") as f:
        prev_query = "NULL"
        best_edit_dist = 99999999
        best_mapq = -1
        best_line = "NULL"
        for line in f:
            curr = line.split("\t")
            if (curr[2] != "*"):
                edit_dist_inc = 11 # tries to find the column that stores edit distance, this column begins with NM
                while (edit_dist_inc < len(curr)):
                    if curr[edit_dist_inc].split(":")[0] == "NM":
                        break
                    edit_dist_inc += 1
                curr_edit_dist = 99999999
                if edit_dist_inc < len(curr):
                    curr_edit_dist = int(curr[edit_dist_inc].split(":")[2])
                if (prev_query == curr[0]): # if we are still on the same read as before, still need to compare
                    if (curr_edit_dist < best_edit_dist): # checking if line has lowest (best) edit dist
                        best_line = line
                        best_edit_dist = curr_edit_dist
                        best_mapq = int(curr[4])
                    elif (curr_edit_dist == best_edit_dist): # checking if line has a higher mapq, if edit dists are equal
                        if (best_mapq != 255):
                            if (int(curr[4]) > best_mapq):
                                best_line = line
                                best_edit_dist = curr_edit_dist # shouldnt matter as they are equal anyways
                                best_mapq = int(curr[4])
                        else: 
                            if (int(curr[4]) != 255): # remember: 255 means unascertainable, hence curr match is best match
                                best_line = line
                                best_edit_dist = curr_edit_dist
                                best_mapq = int(curr[4])
                    prev_query = curr[0]

                else:
                    if (best_line != "NULL"):
                        to_print = best_line.split("\t")
                        read = to_print[0]
                        accession = to_print[2]

                        best_hits_output.write(read + "\t" + accession + "\n")

                    prev_query = curr[0]
                    best_edit_dist = curr_edit_dist
                    best_mapq = int(curr[4])
                    best_line = line

        if (best_line != "NULL"):
            to_print = best_line.split("\t")
            read = to_print[0]
            accession = to_print[2]

            best_hits_output.write(read + "\t" + accession + "\n")

    best_hits_output.close()
