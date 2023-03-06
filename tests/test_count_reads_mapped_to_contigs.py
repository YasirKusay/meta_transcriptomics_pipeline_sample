from meta_transcriptomics_pipeline.finalisation import getContigReadCount

def test_count_reads_mapped_to_contigs():
    pth = "tests/samples/count_reads_mapped_to_contigs/input.txt"

    #read1	cont1
    #read2	cont2
    #read3	cont2
    #read4	cont1
    #read5	cont1
    #read6	cont3
    #read7	cont2

    exp = {}
    exp["cont1"] = 3
    exp["cont2"] = 3
    exp["cont3"] = 1

    actual = getContigReadCount(pth)

    print(actual)

    assert(exp["cont1"] == actual["cont1"])
    assert(exp["cont2"] == actual["cont2"])
    assert(exp["cont3"] == actual["cont3"])

    print("testing finished")