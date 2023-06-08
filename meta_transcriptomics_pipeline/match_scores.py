import os
import pandas as pd
from meta_transcriptomics_pipeline.getScientificNames import getScientificNames
from meta_transcriptomics_pipeline.helpers import run_shell_command

def match_scores(contig_read_mappings, combined_scores, dirpath, outfile, taxids_location):
    temp_out = dirpath + "/temp_out"
    wf = open(temp_out, "w")

    contig_read_mappings_sorted = dirpath + "/contig_read_mappings_sorted"
    run_shell_command("LC_COLLATE=C sort -k1 " + contig_read_mappings + " > " + contig_read_mappings_sorted)
    combined_scores_sorted = dirpath + "/combined_scores_sorted"
    run_shell_command("LC_COLLATE=C sort -k1 " + combined_scores + " > " + combined_scores_sorted)

    all_taxids = []
    with open(contig_read_mappings_sorted, "r") as f:
        for line in f:
            curr = line.split("\t")
            all_taxids.append(curr[2].strip())
    
    sci_names = getScientificNames(all_taxids, taxids_location)

    # this file contains read name \t accession of seq that it mapped onto \t taxid
    mf = open(contig_read_mappings_sorted, "r")
    curr_mf_line = mf.readline().split("\t")
    curr_mf_read = curr_mf_line[0]
    curr_mf_taxid = curr_mf_line[2].strip()

    with open(combined_scores_sorted, "r") as cs:
        for line in cs:
            curr_cs_line = line.split("\t")
            curr_cs_read = curr_cs_line[0]

            if curr_cs_read == curr_mf_read and sci_names[curr_mf_taxid] != "Unknown":
                wf.write(sci_names[curr_mf_taxid].replace(",", " ") + "\t" + curr_cs_line[2] + "\t" + curr_cs_line[3] + "\t" + curr_cs_line[4] + "\t" + curr_cs_line[5].strip() + "\n")

                curr_mf_line = mf.readline().split("\t")
                if (len(curr_mf_line) == 1):
                    break

                curr_mf_read = curr_mf_line[0]
                curr_mf_taxid = curr_mf_line[2].strip()

            elif (curr_mf_read > curr_cs_read):
                continue
            else:
                curr_mf_line = mf.readline().split("\t")
                if (len(curr_mf_line) == 1):
                    break

                curr_mf_read = curr_mf_line[0]
                curr_mf_taxid = curr_mf_line[2].strip()

    wf.close()
    mf.close()

    df = pd.read_csv(temp_out, sep='\t', names=["Taxid", "E-val", "Bitscore", "Percent-id", "Length"])
    df = df.astype({"Taxid": str, "E-val": float, "Bitscore": float, "Percent-id": float, "Length": int})
    df = df.groupby("Taxid").mean().reset_index()
    df.to_csv(outfile, sep="\t", index=False, header=False)

    os.remove(temp_out)