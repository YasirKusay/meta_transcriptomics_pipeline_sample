def filter_result(infile, pid, e_val, bitscore):
    filter_out = []
    with open(infile, "r") as f:
        for line in f:
            curr = line.split()
            if (float(curr[1]) < pid):
                filter_out.append(curr[0])
                continue
            if (float(curr[2]) < e_val):
                filter_out.append(curr[0])
                continue
            if (float(curr[3].strip()) < e_val):
                filter_out.append(curr[0])
                continue 

    return filter_out

def get_filtered_taxids(taxids_to_remove, counts_file, outfile):
    wf = open(outfile, "w")
    with open(counts_file, "r") as f:
        unknown_amount = 0
        for line in f:
            curr = line.strip().split()
            percent = curr[1]
            taxid = curr[0]
            if taxid in taxids_to_remove:
                unknown_amount += float(percent)
            elif taxid == "Unknown":
                unknown_amount += float(percent)
            else:
                wf.write(str(taxid) + "\t" + str(percent) + "\n")
 
        wf.write("Unknown\t" + str(unknown_amount) + "\n")
        
    wf.close()
