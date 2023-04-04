import filecmp
import os
from meta_transcriptomics_pipeline.count_num_lines import countNumSeqs, countNumLines

def test_count_lines():

    # these tests are more of a 'sanity test', i.e. the functions are very
    # simple, however, sometimes for an unexplained reason cythonized 
    # modules throw segmentation faults, which is what we want to
    # guard against

    inp = "tests/samples/test_count_lines/"

    assert(countNumLines(inp + "countLines1.txt") == 4)
    assert(countNumLines(inp + "countLines2.txt") == 3)
    assert(countNumSeqs(inp + "countFasta.fa", True) == 4)
    assert(countNumSeqs(inp + "countFastq.fq", False) == 5)