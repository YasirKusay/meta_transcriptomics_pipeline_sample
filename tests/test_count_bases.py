import filecmp
import os
from meta_transcriptomics_pipeline.main_pipeline import getReadsLength

def test_count_bases():
    inp = "tests/samples/test_count_bases/input.txt"
    exp = "tests/samples/test_count_bases/exp.txt"
    out = "tests/samples/test_count_bases/actual.txt"

    if os.path.exists(out):
        os.remove(out)

    getReadsLength(inp, out, None, True)
    assert(filecmp.cmp(exp, out, shallow=False))

