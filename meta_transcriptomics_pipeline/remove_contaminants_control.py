import subprocess
import os
import math
import re
from helpers import check_fail

def run_kraken(combined, kraken_index, threads, path):
    kraken_outs = []
    i = 0

    for item in combined:
        if os.path.isfile(item) == False:
            print(item + "Does not exist")
            return None

    while i < len(combined):
        # paired end check 
        file1 = combined[i]
        file2 = combined[i+1]

        curr_path = file1.split(".")[0]
        outfile = curr_path.split("/")[-1]
        kraken_outs.append(path + "/" + outfile)

        #kraken_command = "kraken2 --db " + kraken_index + " --threads " + str(threads) + " --output " + outfile + " --paired " + file1 + " " + file2
        #new_command = subprocess.run(kraken_command, shell=True)
        #if check_fail("kraken", new_command, []) is True: return None
        i += 2

    print(kraken_outs)

    return kraken_outs

def run_rcf(taxdump_location, kraken_controls, kraken_outs, num_controls, rcf_out, path):
    rcf_command = "rcf -n " + taxdump_location + " -k " + " -k ".join(kraken_controls) + " -k " +  " -k ".join(kraken_outs) + " -c " + str(len(kraken_controls)) + " -o " + path + "/rcf_out_temp" + " -s KRAKEN -y 25 -d > " + rcf_out
    print(rcf_command)
    new_command = subprocess.run(rcf_command, shell=True)
    if check_fail("kraken", new_command, []) is True: return None
    return True

# controls and sequences are lists
def remove_contaminants_control(controls, sequences, kraken_index, rcf_out, taxdump_location, threads, path):
    # need to firstly run kraken on all 
    kraken_outs = []
    # assume that it is paired end for now
    combined = controls + sequences
 
    print("Running kraken")
    #kraken_outs = run_kraken(combined, kraken_index, threads, path)
    #if kraken_outs is None:
    #    return None

    # now we can run the actual recentrifuge command
    print("Running RCF")
    if run_rcf(taxdump_location, controls, sequences, len(controls), rcf_out, path) is None:
        return None

    # now need to look through file to return the contaminant samples
    # contaminants are stored at the taxid level
    contaminants = []
    print("CONT")
    #contaminant_levels = ["critical", "severe", "mild cont", "just-ctrl", "other cont"]
    # contaminant_levels = ["critical", "severe", "mild cont", "crossover", "other cont"]
    contaminant_levels = ["critical", "severe", "crossover"]
    contaminant_levels_2 = ["mild", "other"]
    with open(rcf_out, "r") as r:
        for line in r:
            for level in contaminant_levels:
                #if (level in line):
                #if re.search("^\^\[\[[0-9]+m" + level, line) is not None:
                if level in line:
                    splits = line.split()
                    contaminants.append(splits[2])
                    print(splits[2])
                    break

            for level in contaminant_levels_2:
                #if (level in line):
                #if re.search("^\^\[\[[0-9]+m" + level + " cont", line) is not None:
                if level + " cont" in line:
                    splits = line.split()
                    contaminants.append(splits[3])
                    print(splits[3])
                    break

    return contaminants 
