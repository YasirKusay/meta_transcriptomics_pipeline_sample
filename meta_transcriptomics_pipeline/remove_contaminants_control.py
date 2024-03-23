import logging
from meta_transcriptomics_pipeline.helpers import run_shell_command

log = logging.getLogger(__name__)

def run_rcf(taxdump_location, kraken_controls_out, kraken_others_out, rcf_out, tmp_path):
    rcf_command = "rcf -n " + taxdump_location + " -k " + " -k ".join(kraken_controls_out) + " -k " +  " -k ".join(kraken_others_out) + " -c " + str(len(kraken_controls_out)) + " -o " + tmp_path + "/rcf_out_temp" + " -s KRAKEN -y 25 -d | sed -e 's/\x1b\[[0-9;]*m//g' > " + rcf_out
    run_shell_command(rcf_command)
    return True

# makes it easier to write tests for it
def identify_contaminants(rcf_outfile):
    contaminants = []
    contaminant_levels = ["critical", "severe", "crossover", "mild cont", "other cont"]
    with open(rcf_outfile, "r") as r:
        for line in r:
            for level in contaminant_levels:
                if level in line:
                    splits = line.split("\t")
                    contaminant_info = splits[1].split(" ")
                    contaminant_taxid = contaminant_info[1].strip()
                    contaminants.append(contaminant_taxid)
                    break

    return contaminants

# controls and others are lists
def remove_contaminants_control(kraken_controls_out, kraken_others_out, rcf_out, taxdump_location, tmp_path):
    # now we can run the actual recentrifuge command
    log.warning("Running RCF")
    if run_rcf(taxdump_location, kraken_controls_out, kraken_others_out, rcf_out, tmp_path) is None:
        return None

    # now need to look through file to return the contaminant samples
    # contaminants are stored at the taxid level
    
    return identify_contaminants(rcf_out) 
