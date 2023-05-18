import os
import operator

def get_abundance(joined, total_reads, n50, output_file, contaminants, bad_taxids):
    taxids = {}
    lengths = {}
    counts = {}

    denom = 0.0

    with open(joined, "r") as f:
        for line in f:
            curr = line.split()
            if (curr[3].strip() in contaminants or curr[3].strip() in bad_taxids):
                total_reads = total_reads - int(curr[2])

    denom = total_reads

    with open(joined, "r") as f:
        for line in f:
            curr = line.split()
            if (int(curr[1]) < 50):
                continue
            if (curr[3].strip() in contaminants or curr[3].strip() in bad_taxids): # skipping contaminants
                continue
            name = curr[0]
            length = int(curr[1])
            if (name.startswith("k") is False): # is not a contig
                length = n50
            count = int(curr[2])
            taxid = int(curr[3])

            taxids[name] = taxid
            lengths[name] = length
            counts[name] = count

    # solution below adapted from RSEM
    # https://github.com/deweylab/RSEM/blob/e4dda70e90fb5eb9b831306f1c381f8bbf71ef0e/WriteResults.h
    # function adapted: calcExpressionValues
    # only difference is I did not use expected length, using the real length instead

    for key in lengths:
        counts[key] = counts[key]/denom

    fpkm = {}
    # calculate FPKM
    for key in lengths:
        fpkm[key] = (counts[key] * 1e9)/lengths[key]

    denom = 0.0
    for key in lengths:
        denom += fpkm[key]

    tpm = {}
    for key in lengths:
        tpm[key] = fpkm[key]/(denom * 1e6)

    for key in tpm:
        print(str(key) + "\t" + str(tpm[key]) + "\n")

    final_tpm = {}
    # calculating tpm but for the taxids
    for key in lengths:
        if taxids[key] in final_tpm.keys():
            final_tpm[taxids[key]] += tpm[key]
        else:
            final_tpm[taxids[key]] = tpm[key]

    if os.path.isfile(output_file):
        os.remove(output_file)
    wf = open(output_file, "w")
    final_tpm_sorted = dict( sorted(final_tpm.items(), key=operator.itemgetter(1),reverse=True))
    for key in final_tpm_sorted:
        if (key not in contaminants or key not in bad_taxids):
            wf.write(str(key) + "\t" + str(final_tpm_sorted[key]) + "\n")
    wf.close()