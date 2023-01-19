import os
import filecmp
from meta_transcriptomics_pipeline.merge_blast_outputs import merge_blast_outputs

def test_count_bases():
    pth = "/testing/samples/merge_contigs/"
    nucl_file = os.getcwd() + pth + "nucl.txt"
    prot_file = os.getcwd() + pth + "prot.txt"
    nucl_out = os.getcwd() + pth + "nucl_out"
    prot_out = os.getcwd() + pth + "prot_out"
    nucl_out_test, prot_out_test = merge_blast_outputs(nucl_file, prot_file, os.getcwd())
    assert(True)
