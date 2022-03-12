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
        stop = False
        for line in cs:
            curr = line.split("\t")
            curr_read = curr[0]
            curr_accession = curr[1]

            if (curr_read == other_read):
                wf.write(other_taxid.strip() + "\t" + curr[2] + "\t" + curr[3] + "\t" + curr[4].strip() + "\n")

                other_line = mf.readline().split("\t")
                if (len(other_line) == 1):
                    stop = True
                    break

                other_read = other_line[0]
                other_accession = other_line[1]
                other_taxid = other_line[2]

            elif (other_read > curr_read):
                continue
                #while (other_read > curr_read):
                #    print("skipping" + other_read)
                #    other_line = mf.readline().split("\t")
                #    print(other_line)
                #    if (len(other_line) == 1):
                #        stop = True
                #        break
                #    other_read = other_line[0]
                #    other_accession = other_line[1]
                #    other_taxid = other_line[2]
     
                #if (stop is True):
                #    break

                #if (curr_read == curr_accession):
                #    wf.write(other_taxid.strip() + "\t" + curr[2] + "\t" + curr[3] + "\t" + curr[4].strip() + "\n")

    wf.close()
    mf.close()

    temp_out_sorted = dirpath + "/temp_out_sorted"
    new_command = subprocess.run("LC_COLLATE=C sort -k1 " + temp_out + " > " + temp_out_sorted, shell=True)

    df = pd.read_csv(temp_out, sep='\t', names=["Taxid", "E-val", "Bitscore", "Pid"])
    df = df.astype({"Taxid": str, "Pid": float, "E-val": float, "Bitscore": float})
    df = df.groupby("Taxid").mean().reset_index()
    df.to_csv(outfile, sep="\t", index=False, header=False)
