import filecmp
from meta_transcriptomics_pipeline.main_pipeline import getReadsLength

def test_count_bases():
    inp = "tests/samples/test_count_bases/input.txt"
    exp = ""
    out = "tests/samples/test_count_bases/actual.txt"
    getReadsLength(inp, out, None, True)
    assert(filecmp.cmp(out, ))
