import subprocess

def adjust_tpms(read1, read2, tpm_file):
    cmd = "cat " +  tpm_file + " | awk 'NR > 1' | sort -n -k 6 | uniq"
    