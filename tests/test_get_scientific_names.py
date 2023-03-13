import os

from meta_transcriptomics_pipeline.getScientificNames import getScientificNames

def test_get_scientific_names():
    names="tests/samples/get_scientific_names/"

    exp1 = {"1": "Species A", "2": "Species B", "5": "Species E", "7": "Species G", "8": "8", "9": "9", "10": "Species J", "Unknown": "Unknown"}
    assert(exp1 == getScientificNames(["1", "2", "5", "7", "8", "9", "10"], names))

    exp2 = {"1": "Species A", "2": "Species B", "5": "Species E", "7": "Species G", "8": "8", "9": "9", "10": "Species J", "11": "11", "12": "12", "Unknown": "Unknown"}
    assert(exp2 == getScientificNames(["1", "2", "5", "7", "8", "9", "10", "11", "12"], names))

    exp2 = {"1": "Species A", "2": "Species B", "Unknown": "Unknown"}
    assert(exp2 == getScientificNames(["1", "2"], names))

    exp2 = {"1": "Species A", "2": "Species B", "8": "8", "Unknown": "Unknown"}
    assert(exp2 == getScientificNames(["1", "2", "8"], names))