import subprocess
import pandas as pd
from ete3 import NCBITaxa
from meta_transcriptomics_pipeline.getScientificNames import getScientificNames

def match_scores(mappings, combined_scores, dirpath, outfile, taxids_location):
    temp_out = dirpath + "/temp_out"
    wf = open(temp_out, "w")

    mappings_sorted = dirpath + "/mappings_sorted"
    new_command = subprocess.run("LC_COLLATE=C sort -k1 " + mappings + " > " + mappings_sorted, shell=True)
    combined_scores_sorted = dirpath + "/combined_scores_sorted"
    new_command = subprocess.run("LC_COLLATE=C sort -k1 " + combined_scores + " > " + combined_scores_sorted, shell=True)

    mf = open(mappings_sorted, "r")
    other_line = mf.readline().split("\t")
    other_read = other_line[0]
    other_accession = other_line[1]
    other_taxid = other_line[2]

    ncbi = NCBITaxa()

    with open(combined_scores_sorted, "r") as cs:
        stop = False
        for line in cs:
            curr = line.split("\t")
            curr_read = curr[0]
            curr_accession = curr[1]

            if (curr_read == other_read):

                lineage = None
                try:
                    lineage = ncbi.get_lineage(other_taxid.strip())
                except:
                    print(str(other_taxid.strip()) + " NOT FOUND")
                    continue

                lineage2ranks = ncbi.get_rank(lineage)
                ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
                r = ranks2lineage.get("species", "Unknown")
                wf.write(str(r) + "\t" + curr[2] + "\t" + curr[3] + "\t" + curr[4].strip() + "\n")

                other_line = mf.readline().split("\t")
                if (len(other_line) == 1):
                    stop = True
                    break

                other_read = other_line[0]
                other_accession = other_line[1]
                other_taxid = other_line[2]

            elif (other_read > curr_read):
                continue
            else:
                other_line = mf.readline().split("\t")
                if (len(other_line) == 1):
                    stop = True
                    break

                other_read = other_line[0]
                other_accession = other_line[1]
                other_taxid = other_line[2]

    wf.close()
    mf.close()

    temp_out_sorted = dirpath + "/temp_out_sorted"
    new_command = subprocess.run("LC_COLLATE=C sort -k1 " + temp_out + " > " + temp_out_sorted, shell=True)

    df = pd.read_csv(temp_out, sep='\t', names=["Taxid", "E-val", "Bitscore", "Pid"])
    df = df.astype({"Taxid": str, "Pid": float, "E-val": float, "Bitscore": float})
    df = df.groupby("Taxid").mean().reset_index()
    temp_out_2 = dirpath + "/temp_out_2"
    df.to_csv(temp_out_2, sep="\t", index=False, header=False)

    taxids = []
    with open(temp_out_2, "r") as mf:
        for line in mf:
            curr = line.split("\t")
            if (curr[0] != "Unknown"):
                taxids.append(int(curr[0]))
   
    sci_names = getScientificNames(taxids, taxids_location)

    fout = open(outfile, "w")

    with open(temp_out_2, "r") as mf:
        for line in mf:
            curr = line.split("\t")
            fout.write(sci_names[curr[0]] + "\t" + curr[1] + "\t" + curr[2] + "\t" + curr[3].strip() + "\n")
