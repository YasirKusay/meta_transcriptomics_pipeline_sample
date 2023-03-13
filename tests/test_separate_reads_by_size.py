import filecmp
import os

from meta_transcriptomics_pipeline.separate_reads_by_size import separate_reads_by_size

def test_separate_reads_by_size():
    pth="tests/samples/separate_reads_by_size/"

    if os.path.exists(pth + "outLong.fq"):
        os.remove(pth + "outLong.fq")
        
    if os.path.exists(pth + "outShort.fq"):
        os.remove(pth + "outShort.fq")

    separate_reads_by_size(pth + "inputReads.fq", pth + "outLong.fq", pth + "outShort.fq")

    assert(filecmp.cmp(pth + "expShort.fq", pth + "outShort.fq", shallow=False) == True)
    assert(filecmp.cmp(pth + "expLong.fq", pth + "outLong.fq", shallow=False) == True)