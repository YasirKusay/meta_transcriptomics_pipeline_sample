import os
import filecmp
from meta_transcriptomics_pipeline.merge_blast_outputs import merge_blast_outputs

def test_count_bases():
    pth = "tests/samples/merge_blast_outputs"
    nucl_file = pth + "/nucl.txt"
    prot_file = pth + "/prot.txt"

    nucl_exp = pth + "/nucl_exp.txt"
    prot_exp = pth + "/prot_exp.txt"

    nucl_out, prot_out = merge_blast_outputs(nucl_file, prot_file, pth)

    assert(filecmp.cmp(nucl_exp, nucl_out, shallow=False))
    assert(filecmp.cmp(prot_exp, prot_out, shallow=False))

    os.remove(pth + "/temp_sam")
    os.remove(pth + "/temp_diamond")
    os.remove(pth + "/snap_diamond_combined_file")
    os.remove(pth + "/snap_diamond_sorted_file")