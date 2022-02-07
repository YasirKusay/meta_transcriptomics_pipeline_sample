import subprocess
from helpers import check_fail

def run_kraken(combined, kraken_index, threads):
    kraken_outs = []
    i = 0
    while i < len(combined):
        # paired end check 
        file1 = combined[i]
        file2 = combined[i+1]

        outfile = file1.split(".")[0]
        kraken_outs.append(outfile)

        kraken_command = "kraken --db " + kraken_index + " --threads " + str(threads) + " --output " + outfile + " --paired " + file1 + " " + file2
        new_command = subprocess.run(kraken_command, shell=True)
        if check_fail("kraken", new_command, []) is True: return None
        i += 2

    return kraken_outs

def run_rcf(taxdump_location, kraken_outs, num_controls, rcf_out):
    rcf_command = "rcf -n " + taxdump_location + " -k " +  " -k ".join(kraken_outs) + " -c " + str(num_controls) + " " + "rcf_out_temp" + " -s KRAKEN -y 25 -d > " + rcf_out
    new_command = subprocess.run(rcf_command, shell=True)
    if check_fail("kraken", new_command, []) is True: return None
    return True

# controls and sequences are lists
def remove_contaminants_control(controls, sequences, kraken_index, rcf_out, taxdump_location, threads):
    # need to firstly run kraken on all 
    kraken_outs = []
    # assume that it is paired end for now
    combined = controls + sequences

    kraken_outs = run_kraken(combined, kraken_index, threads)
    if kraken_outs is None:
        return None

    # now we can run the actual recentrifuge command
    if run_rcf(taxdump_location, kraken_outs, len(controls), rcf_out) is None:
        return None

    # now need to look through file to return the contaminant samples
    # contaminants are stored at the taxid level
    contaminants = []
    #contaminant_levels = ["critical", "severe", "mild cont", "just-ctrl", "other cont"]
    contaminant_levels = ["critical", "severe", "mild cont", "crossover", "other cont"]
    with open(rcf_out, "r") as r:
        for line in r:
            for level in contaminant_levels:
                if (level in line):
                    splits = line.split()
                    other = splits[1].split(" ")
                    contaminants.append(other[0])
                    break

    return contaminants 