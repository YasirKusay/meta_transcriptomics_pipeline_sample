'''
MIT License

Copyright (c) 2018 Chan Zuckerberg Initiative, LLC.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

# originally from https://github.com/chanzuckerberg/czid-workflows/blob/72fa52585b318f29eea24c91496df6e116262315/short-read-mngs/idseq_utils/idseq_utils/paf2blast6.py

import pandas as pd
import os.path
import math
import re
import numpy
from meta_transcriptomics_pipeline.helpers import run_shell_command

nonmatch_pattern = re.compile(r"NM:i:(\d+);")
cigar_pattern = re.compile(r"cg:Z:([A-Za-z0-9]+)")


class QualityCalculations:
    # THE GENOME SIZE IS HARDCODED AND WILL NEED TO BE UPDATED EVERY SINGLE TIME WE USE A NEW INDEX
    # THIS IS THE CURRENT INDEX SIZE AS OF MAY 2023
    def __init__(self, genome_size=1138726549535): 
        self._k = 0.1
        self._lambda = 1.58
        self.genome_size = genome_size

    def calc_bitscore(self, alen, nonmatch):
        score = alen - 2 * nonmatch
        return (score * self._lambda - math.log(self._k)) / math.log(2.0)

    def calc_evalue(self, alen, nonmatch):
        score = alen - 2 * nonmatch
        return self._k * alen * self.genome_size * numpy.exp(-self._lambda * score)

    def calc_gap_openings(self, cigar):
        go = 0
        for char in cigar:
            if char == "I" or char == "D":
                go += 1
        return go


def standardize_paf(paf_file):
    """paf files can have variable number of optional key-value pairs"""
    run_shell_command("sed 's/\t/;/13g' -i {}".format(paf_file))


def main(paf_file, dirpath):
    name = os.path.basename(paf_file.replace(".paf", ""))
    standardize_paf(paf_file)
    df = pd.read_csv(
        paf_file,
        delimiter="\t",
        header=None,
        names=[
            "qname",
            "qlen",
            "qstart",
            "qend",
            "strand",
            "tname",
            "tlen",
            "tstart",
            "tend",
            "nmatch",
            "alen",
            "mapq",
            "other",
        ],
    )
    df["nonmatch"] = df.other.map(
        lambda x: int(re.search(nonmatch_pattern, x).group(1))
    )
    qc = QualityCalculations()
    df["gap_openings"] = df.other.map(
        lambda x: qc.calc_gap_openings(re.search(cigar_pattern, x).group(1))
    )
    df["bitscore"] = [
        qc.calc_bitscore(a, n) for a, n in zip(df["alen"], df["nonmatch"])
    ]
    df["evalue"] = [qc.calc_evalue(a, n) for a, n in zip(df["alen"], df["nonmatch"])]
    df["percent_ident"] = [
        (nmatch / a) * 100 for a, nmatch in zip(df["alen"], df["nmatch"])
    ]
    df = df.round({"bitscore": 3, "percent_ident": 3})
    blast = df.loc[
        :,
        [
            "qname",
            "tname",
            "percent_ident",
            "alen",
            "nonmatch",
            "gap_openings",
            "qstart",
            "qend",
            "tstart",
            "tend",
            "evalue",
            "bitscore",
            "qlen"
        ],
    ]

    blast["qstart"] = blast["qstart"] + 1
    blast["tstart"] = blast["tstart"] + 1

    blast.to_csv(f"{dirpath}/{name}_frompaf.m8", sep="\t", index=None, header=False)


def paf2blast6(paf_file, dirpath):
    main(paf_file, dirpath)
