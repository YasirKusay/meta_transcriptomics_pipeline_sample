# file takes the reads mapped to contigs sam file,
# retrieves the contig that is mapped to their read
# further retrieves the taxid they are mapped to

def map_reads_to_contigs(reads_mapped_to_contigs, outfile):
    to_write = open(outfile, "w")

    with open(reads_mapped_to_contigs, "r") as r:
        for line in r:
            curr = line.split()

            if (len(curr) <= 3):
                next

            if (curr[2] != "*"):
                read = curr[0]
                contig = curr[2]

                to_write.write(read + "\t" + contig + "\n")
