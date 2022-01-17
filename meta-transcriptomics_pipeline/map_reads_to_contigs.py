# file takes the reads mapped to contigs sam file,
# retrieves the contig that is mapped to their read
# further retrieves the taxid they are mapped to
import os
import re

def map_reads_to_contigs(reads_mapped_to_contigs, outfile):
    
    if os.path.isfile(outfile):
        print("DELETED")
        os.remove(outfile)
    to_write = open(outfile, "w")

    with open(reads_mapped_to_contigs, "r") as r:
        already_seen = [] # temp solution, delete later
        num = 1 # also temp solution
        for line in r:
            curr = line.split()

            #if (len(curr) <= 3):
            #    continue

            if (len(re.findall("^@", curr[0])) > 0):
                continue

            if num == 2:
                num = 1
                continue
            num += 1

            if (curr[2] != "*"):
                read = curr[0]
                contig = curr[2]
                #print(read)
                #if read in already_seen:
                #    continue
                #else:
                #    already_seen.append(read)
                to_write.write(read + "\t" + contig + "\n")

    to_write.close()
    with open(outfile, "r") as r:
        num = 1 # also temp solution
        count = 0
        for line in r:
            count += 1

        print("printing count: " + str(count))
