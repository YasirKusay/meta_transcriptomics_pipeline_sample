import os
import filecmp

from meta_transcriptomics_pipeline.merge_sams import merge_sams

def test_count_bases():
    nucl_file = "tests/samples/merge_contigs/nucl.sam"
    prot_file = "tests/samples/merge_contigs/prot.sam"
    nucl_out = "tests/samples/merge_contigs/nucl_out.txt"
    prot_out = "tests/samples/merge_contigs/prot_out.txt"
    nucl_exp = "tests/samples/merge_contigs/nucl_exp.txt"

    if os.path.exists(nucl_out):
        os.remove(nucl_out)
    if os.path.exists(prot_out):
        os.remove(prot_out)

    print(os.getcwd())

    merge_sams(nucl_file, prot_file, os.getcwd(), nucl_out, prot_out)
    assert(filecmp.cmp(nucl_exp, nucl_out, shallow=False))

    if os.path.exists("snap_diamond_combined_file"):
        os.remove("snap_diamond_combined_file")
    if os.path.exists("snap_diamond_sorted_file"):
        os.remove("snap_diamond_sorted_file")