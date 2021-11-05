# file takes the 2 sam files, one from snap and the other for diamond and merges them
# and finds the best hit possible
# best hits are decided by edit distance (1) and mapq (2)

import sys
import os
import subprocess
from helpers import check_fail

def join(file, mapping_file, join_file):
    command = "join -1 2 -2 1 -o \"1.1, 1.2, 2.1\" <(sort -k2 " + file + ") <(sort -k1 " + mapping_file + ") >> " + join_file
    new_command = subprocess.run(command, shell=True)
    if check_fail("join", new_command, []) is True: return None
    
def join_taxid_contigs(output_snap, output_diamond, nucl_map, prot_map, path):
    final_out = path + "/nucl_prot_taxids.txt"
    if join(output_snap, nucl_map, final_out) == None or join(output_diamond, prot_map, final_out) == None: 
        return None

    return final_out
    