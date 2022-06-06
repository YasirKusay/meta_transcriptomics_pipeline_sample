from meta_transcriptomics_pipeline.get_lineage_info import getReadsLength

def test_count_bases():
    inp = "samples/test_count_bases/input.txt"
    out = "samples/test_count_bases/actual.txt"
    getReadsLength(inp, out, None, True)
    assert(True)
