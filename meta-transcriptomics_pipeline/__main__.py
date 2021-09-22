import argparse
import sys 
import multiprocessing
from main_pipeline import run_pipeline

def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect contamination in an assembly"
    )

    parser.add_argument(
        "inp1",
        type=str,
        help="read1 input name"
    )

    parser.add_argument(
        "inp2",
        type=str,
        help="read1 input name"
    )

    parser.add_argument(
        "sortmerna_index",
        type=str,
        help="sortmerna-index"
    )

    parser.add_argument(
        "--qualified_quality_phred",
        type=int,
        default=15,
        help=" DURING QC. The quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified"
    )

    parser.add_argument(
        "--unqualified_percent_limit",
        type=int,
        default=45,
        help=" DURING QC. How many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%"
    )

    parser.add_argument(
        "--average_qual",
        type=int,
        default=0,
        help=" DURING QC. if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])"
    )

    parser.add_argument(
        "--length_required",
        type=int,
        default=15,
        help="DURING QC. Reads shorter than length_required will be discarded, default is 15. (int [=15])"
    )

    parser.add_argument(
        "--complexity_threshold",
        type=int,
        default=30,
        help="DURING QC. The threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )
    # should we have options for frameshift, gapopen, gapextend penalties, and scoring matrix
    # also an option for tantan repeat masking, composition matrix

    # will we have access to a taxonmap/taxon names for make db


    parser.set_defaults(func=run_pipeline)
    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
