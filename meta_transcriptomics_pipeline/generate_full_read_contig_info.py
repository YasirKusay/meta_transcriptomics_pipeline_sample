def generate_full_read_contig_info(combined_blast_scores_sorted, read_contigs_mappings, taxid_lineages, outfile):

    wf = open(outfile, "w")

    mf = open(read_contigs_mappings, "r")
    curr_mf_line = mf.readline().strip().split("\t")
    curr_mf_read = curr_mf_line[0]
    curr_mf_taxid = curr_mf_line[2]

    for key in taxid_lineages.keys():
        taxid_lineages[key].reverse()

    # the combined_blast_scores_sorted file should always include the reads in read_contigs_mappings but not the other way around
    # this algorithm will ensure that we only fetch the reads from combined_blast_scores_sorted that are in read_contigs_mappings

    with open(combined_blast_scores_sorted, "r") as cbs:
        for line in cbs:
            curr_blast_line_split = line.strip().split("\t")
            curr_blast_read = curr_blast_line_split[0]

            if curr_blast_read == curr_mf_read and curr_mf_taxid in taxid_lineages.keys():
                # curr_blast_line_split[:-3] if you want to exclude query length/whether the match was a read/protein
                wf.write(line.strip() + "\t" + "\t".join(taxid_lineages[curr_mf_taxid]) + "\n")

                # files in the combined_blast_scores_sorted file may contain duplicate reads
                # immediately fetch the next mf line to prevent the same read being added multiple
                # times as we are forced to skip this line

                curr_mf_line = mf.readline().strip().split("\t")
                if (len(curr_mf_line) == 1):
                    break
                curr_mf_read = curr_mf_line[0]
                curr_mf_taxid = curr_mf_line[2]

            elif (curr_mf_read > curr_blast_read):
                continue
            else:
                curr_mf_line = mf.readline().strip().split("\t")
                if (len(curr_mf_line) == 1):
                    break

                curr_mf_read = curr_mf_line[0]
                curr_mf_taxid = curr_mf_line[2]

    wf.close()