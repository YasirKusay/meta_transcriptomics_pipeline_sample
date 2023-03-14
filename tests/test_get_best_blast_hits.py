import os
import filecmp
from meta_transcriptomics_pipeline.get_best_blast_hits import get_best_blast_hits

def test_count_bases():
    pth = "tests/samples/get_best_blast_hits"
    nucl_file = pth + "/nucl.txt"
    prot_file = pth + "/prot.txt"

    nucl_exp = pth + "/nucl_exp.txt"
    prot_exp = pth + "/prot_exp.txt"

    nucl_out = pth + "/nuclOut.txt"
    prot_out = pth + "/protOut.txt"

    if os.path.exists(pth + "/analysis") is False:
        os.mkdir(pth + "/analysis")

    get_best_blast_hits(nucl_file, prot_file, pth, nucl_out, prot_out)

    assert(filecmp.cmp(nucl_exp, nucl_out, shallow=False))
    assert(filecmp.cmp(prot_exp, prot_out, shallow=False))

    os.remove(pth + "/analysis/nt_nr_alignments_combined")
    os.remove(pth + "/analysis/nt_nr_alignments_combined_sorted")
    os.remove(nucl_out)
    os.remove(prot_out)