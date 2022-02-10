# file takes the reads mapped to contigs sam file,
# retrieves the contig that is mapped to their read
# further retrieves the taxid they are mapped to
import os
from merge_sams import merge_sams

def map_reads_to_contigs(reads_mapped_to_contigs, outfile, path):
    
    if os.path.isfile(outfile):
        os.remove(outfile)
    to_write = open(outfile, "w")

    empty_protein_file = path + "/empty_protein_file"

    open(empty_protein_file, "a").close() # creating an empty list for the function
    protein_out = path + "/protein_out"

    final_contig_read_sam_file = path + "/final_contig_read_sam_file"

    # picking best contig-read assignment based on edit distance and mapping quality
    merge_sams(reads_mapped_to_contigs, empty_protein_file, path, snap_out = final_contig_read_sam_file, diamond_out = protein_out)

    with open(final_contig_read_sam_file, "r") as r:
        for line in r:
            curr = line.split()

            read = curr[0]
            contig = curr[1]
            to_write.write(read + "\t" + contig + "\n")

    to_write.close()
