import argparse
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
        help="path to the sortmerna index"
    )

    parser.add_argument(
        "snap_human_index",
        type=str,
        help="path to the snap human index"
    )

    parser.add_argument(
        "kraken_db",
        type=str,
        help="kraken_db_location"
    )
    
    parser.add_argument(
        "minimap2_index",
        type=str,
        help="path to the minimap2 index to align reads/contigs"
    )

    parser.add_argument(
        "diamond_index",
        type=str,
        help="path to the diamond index to align reads/contigs"
    )

    parser.add_argument(
        "dirpath",
        type=str,
        help="dirpath"
    )

    parser.add_argument(
        "nucl_accession_taxid_mapping",
        type=str,
        help="nucl_accession_taxid_mapping"
    )

    parser.add_argument(
        "prot_accession_taxid_mapping",
        type=str,
        help="prot_accession_taxid_mapping"
    )

    parser.add_argument(
        "taxdump_location",
        type=str,
        help="taxdump_location"
    )

    parser.add_argument(
        "--control_sequences",
        type=str,
        nargs="+",
        help="control sequences for contamination detection"
    )

    parser.add_argument(
        "--other_sequences",
        type=str,
        nargs="+",
        help="normal sequences for contamination detection"
    )

    parser.add_argument(
        "--qualified_quality_phred",
        type=str,
        default="15",
        help=" DURING QC. The quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified"
    )

    parser.add_argument(
        "--unqualified_percent_limit",
        type=str,
        default="45",
        help=" DURING QC. How many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%"
    )

    parser.add_argument(
        "--average_qual",
        type=str,
        default="0",
        help=" DURING QC. if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])"
    )

    parser.add_argument(
        "--length_required",
        type=str,
        default="50",
        help="DURING QC. Reads shorter than length_required will be discarded, default is 15. (int [=15])"
    )

    parser.add_argument(
        "--complexity_threshold",
        type=str,
        default="30",
        help="DURING QC. The threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])"
    ) 

    parser.add_argument(
        "--pid_filter",
        type=float,
        default=0,
        help="pid_filter"
    )

    parser.add_argument(
        "--evalue_filter",
        type=float,
        default=0,
        help="evalue_filter"
    )


    parser.add_argument(
        "--bitscore_filter",
        type=float,
        default=0,
        help="bitscore_filter"
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
