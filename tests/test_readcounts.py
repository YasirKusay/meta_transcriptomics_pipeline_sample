import filecmp
import os

from meta_transcriptomics_pipeline.finalisation import countReads

# to counteract the float point error, just consider 1 decimal place for first test and 2 decimal places for the second test
def fix_decimals(infile, outfile, decimal_points):
    wf = open(outfile, "w")
    with open(infile, "r") as f:
        for line in f:
            curr = line.split("\t")
            curr[1] = float(curr[1].strip())
            curr[1] = str(round(curr[1], decimal_points))
            curr[1] = curr[1] + "\n"
            wf.write("\t".join(curr))

    wf.close()

def test_readcounts():
    if os.path.exists("tests/samples/count_reads/prelim.txt"):
        os.remove("tests/samples/count_reads/prelim.txt")
        
    if os.path.exists("tests/samples/count_reads/prelim2.txt"):
        os.remove("tests/samples/count_reads/prelim2.txt")

    if os.path.exists("tests/samples/count_reads/actual.txt"):
        os.remove("tests/samples/count_reads/actual.txt")
        
    if os.path.exists("tests/samples/count_reads/actual2.txt"):
        os.remove("tests/samples/count_reads/actual2.txt")

    countReads("tests/samples/count_reads/input.txt", 200, "tests/samples/count_reads/prelim.txt", None, [])
    countReads("tests/samples/count_reads/input.txt", 400, "tests/samples/count_reads/prelim2.txt", None, [])

    fix_decimals("tests/samples/count_reads/prelim.txt", "tests/samples/count_reads/actual.txt", 1)
    fix_decimals("tests/samples/count_reads/prelim2.txt", "tests/samples/count_reads/actual2.txt", 2)

    assert(filecmp.cmp("tests/samples/count_reads/actual.txt", "tests/samples/count_reads/expected.txt", shallow=False) == True)
    assert(filecmp.cmp("tests/samples/count_reads/actual2.txt", "tests/samples/count_reads/expected2.txt", shallow=False) == True)