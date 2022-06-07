import argparse
import multiprocessing
from main_pipeline import run_pipeline
from blast_sr import blast_sr
from diamond import diamond
from finalisation import finalisation
from minimap2_contigs import minimap2_contigs
from minimap2_reads import minimap2_reads
from preprocessing import preprocessing

def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect contamination in an assembly"
    )

    subparsers = parser.add_subparsers()

    parser_sequential = subparsers.add_parser(
        "sequential", help="Run the pipeline from the start to the end (sequentially)"
    )

    parser_sequential.add_argument(
        "inp1",
        type=str,
        help="read1 input name"
    )

    parser_sequential.add_argument(
        "inp2",
        type=str,
        help="read1 input name"
    )

    parser_sequential.add_argument(
        "sortmerna_index",
        type=str,
        help="path to the sortmerna index"
    )

    parser_sequential.add_argument(
        "snap_human_index",
        type=str,
        help="path to the snap human index"
    )

    parser_sequential.add_argument(
        "kraken_db",
        type=str,
        help="kraken_db_location"
    )
    
    parser_sequential.add_argument(
        "minimap2_index",
        type=str,
        help="path to the minimap2 index to align reads/contigs"
    )

    parser_sequential.add_argument(
        "diamond_index",
        type=str,
        help="path to the diamond index to align reads/contigs"
    )

    parser_sequential.add_argument(
        "dirpath",
        type=str,
        help="where to place the intermediate outputs"
    )

    parser_sequential.add_argument(
        "nucl_accession_taxid_mapping",
        type=str,
        nargs="+",
        help="nucl_accession_taxid_mapping"
    )

    parser_sequential.add_argument(
        "prot_accession_taxid_mapping",
        type=str,
        nargs="+",
        help="prot_accession_taxid_mapping"
    )

    parser_sequential.add_argument(
        "taxdump_location",
        type=str,
        help="taxdump_location"
    )

    parser_sequential.add_argument(
        "--control_sequences",
        type=str,
        nargs="+",
        help="control sequences for contamination detection"
    )

    parser_sequential.add_argument(
        "--other_sequences",
        type=str,
        nargs="+",
        help="normal sequences for contamination detection"
    )

    parser_sequential.add_argument(
        "--qualified_quality_phred",
        type=str,
        default="15",
        help=" DURING QC. The quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified"
    )

    parser_sequential.add_argument(
        "--unqualified_percent_limit",
        type=str,
        default="45",
        help=" DURING QC. How many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%"
    )

    parser_sequential.add_argument(
        "--average_qual",
        type=str,
        default="0",
        help=" DURING QC. if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])"
    )

    parser_sequential.add_argument(
        "--length_required",
        type=str,
        default="50",
        help="DURING QC. Reads shorter than length_required will be discarded, default is 15. (int [=15])"
    )

    parser_sequential.add_argument(
        "--complexity_threshold",
        type=str,
        default="30",
        help="DURING QC. The threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])"
    ) 

    parser_sequential.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_sequential.set_defaults(func=run_pipeline)

    parser_preprocessing = subparsers.add_parser(
        "preprocessing", help="Perform the preprocessing step"
    )

    ################# PREPROCESSING ###################

    parser_preprocessing.add_argument(
        "inp1",
        type=str,
        help="read1 input name"
    )

    parser_preprocessing.add_argument(
        "inp2",
        type=str,
        help="read1 input name"
    )

    parser_preprocessing.add_argument(
        "sortmerna_index",
        type=str,
        help="path to the sortmerna index"
    )

    parser_preprocessing.add_argument(
        "snap_human_index",
        type=str,
        help="path to the snap human index"
    )

    parser_preprocessing.add_argument(
        "kraken_db",
        type=str,
        help="kraken_db_location"
    )

    parser_preprocessing.add_argument(
        "dirpath",
        type=str,
        help="where to place the intermediate outputs"
    )

    parser_preprocessing.add_argument(
        "--qualified_quality_phred",
        type=str,
        default="15",
        help=" DURING QC. The quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified"
    )

    parser_preprocessing.add_argument(
        "--unqualified_percent_limit",
        type=str,
        default="45",
        help=" DURING QC. How many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%"
    )

    parser_preprocessing.add_argument(
        "--average_qual",
        type=str,
        default="0",
        help=" DURING QC. if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])"
    )

    parser_preprocessing.add_argument(
        "--length_required",
        type=str,
        default="50",
        help="DURING QC. Reads shorter than length_required will be discarded, default is 15. (int [=15])"
    )

    parser_preprocessing.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_preprocessing.set_defaults(func=preprocessing)

    parser_minimap_contigs = subparsers.add_parser(
        "minimap_contig_alignment", help="Perform the minimap contig alignment"
    )

    parser_minimap_contigs.add_argument(
        "minimap2_index",
        type=str,
        help="path to the minimap2 index to align reads/contigs"
    )

    parser_minimap_contigs.add_argument(
        "dirpath",
        type=str,
        help="where to place the intermediate outputs"
    )

    parser_minimap_contigs.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_minimap_contigs.set_defaults(func=minimap2_contigs)

    parser_minimap_reads = subparsers.add_parser(
        "minimap_long_read_alignment", help="Perform the minimap longer read alignment"
    )

    parser_minimap_reads.add_argument(
        "minimap2_index",
        type=str,
        help="path to the minimap2 index to align reads/contigs"
    )

    parser_minimap_reads.add_argument(
        "dirpath",
        type=str,
        help="where to place the intermediate outputs"
    )

    parser_minimap_reads.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_minimap_reads.set_defaults(func=minimap2_reads)

    parser_blast_reads = subparsers.add_parser(
        "blast_short_read_alignment", help="Perform the blast shorter read alignment"
    )

    parser_blast_reads.add_argument(
        "dirpath",
        type=str,
        help="where to place the intermediate outputs"
    )

    parser_blast_reads.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_blast_reads.set_defaults(func=blast_sr)

    parser_diamond = subparsers.add_parser(
        "diamond_alignment", help="Perform the diamond alignment against all sequences"
    )

    parser_diamond.add_argument(
        "dirpath",
        type=str,
        help="where to place the intermediate outputs"
    )

    parser_diamond.add_argument(
        "diamond_index",
        type=str,
        help="path to the diamond index to align reads/contigs"
    )

    parser_diamond.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_diamond.set_defaults(func=diamond)

    parser_finalisation = subparsers.add_parser(
        "finalisation", help="Perform taxnonomic identification, decontamination (optional) and quantification"
    )

    parser_finalisation.add_argument(
        "dirpath",
        type=str,
        help="where to place the intermediate outputs"
    )

    parser_finalisation.add_argument(
        "nucl_accession_taxid_mapping",
        type=str,
        nargs="+",
        help="nucl_accession_taxid_mapping"
    )

    parser_finalisation.add_argument(
        "prot_accession_taxid_mapping",
        type=str,
        nargs="+",
        help="prot_accession_taxid_mapping"
    )

    parser_finalisation.add_argument(
        "taxdump_location",
        type=str,
        help="taxdump_location"
    )

    parser_finalisation.add_argument(
        "--kraken_db",
        type=str,
        help="kraken_db_location"
    )

    parser_finalisation.add_argument(
        "--control_sequences",
        type=str,
        nargs="+",
        help="control sequences for contamination detection"
    )

    parser_finalisation.add_argument(
        "--other_sequences",
        type=str,
        nargs="+",
        help="normal sequences for contamination detection"
    )

    parser_finalisation.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser_finalisation.set_defaults(func=finalisation)
    
    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
