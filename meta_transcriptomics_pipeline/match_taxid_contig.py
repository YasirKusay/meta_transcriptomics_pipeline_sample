import subprocess
from helpers import check_fail

def merge_taxid_out(mapping_file, output_file, join_file):
    command = "join -1 1 -2 1 -o \"2.1, 1.2, 2.3, 2.4, 2.5, 2.6, 2.7\" <(sort -k1 " + mapping_file + ") <(sort -k1 " + output_file + ") > " + join_file
    new_command = subprocess.run(command, shell=True)
    if check_fail("join", new_command, []) is True: return None


    cat out.genes.results | awk 'NR > 1' | sort -n -k 6 | uniq | head -n 1000