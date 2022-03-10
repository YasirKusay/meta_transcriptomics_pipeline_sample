import subprocess
import pandas as pd

def match_scores(mappings, combined_scores, dirpath, outfile):
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

    with open(combined_scores_sorted, "r") as cs:
        for line in cs:
            curr = line.split("\t")
            curr_read = curr[0]
            curr_accession = curr[1]

            if (curr_read == curr_accession):
                wf.write(other_taxid + "\t" + curr[2] + "\t" + curr[10] + "\t" + curr[11] + "\n")

            elif (other_read > curr_accession):
                while (other_read > curr_accession):
                    other_line = mf.readline().split("\t")
                    other_read = other_line[0]
                    other_accession = other_line[1]
                    other_taxid = other_line[2]

                if (curr_read == curr_accession):
                    wf.write(other_taxid + "\t" + curr[2] + "\t" + curr[10] + "\t" + curr[11] + "\n")

    df = pd.read_csv(temp_out, sep='\t', names=["Taxid", "Pid", "E-val", "Bitscore"])
    df.groupby("Taxid").mean().reset_index()
    df.to_csv(outfile, sep="\t", index=False, header=False)
