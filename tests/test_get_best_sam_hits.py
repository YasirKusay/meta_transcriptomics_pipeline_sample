import os
import filecmp

from meta_transcriptomics_pipeline.get_best_sam_hits import get_best_sam_hits

def test_count_bases():
    nucl_file = "tests/samples/get_best_sam_hits/nucl.sam"
    nucl_out = "tests/samples/get_best_sam_hits/nucl_out.txt"
    nucl_exp = "tests/samples/get_best_sam_hits/nucl_exp.txt"

    if os.path.exists(nucl_out):
        os.remove(nucl_out)

    get_best_sam_hits(nucl_file, nucl_out,  os.getcwd())
    assert(filecmp.cmp(nucl_exp, nucl_out, shallow=False))

    if os.path.exists("snap_sorted_file"):
        os.remove("snap_sorted_file")