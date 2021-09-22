import argparse
import sys 
import multiprocessing
from main_pipeline import detection

def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect contamination in an assembly"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help="the number of cpu threads to use"
    )

    parser.add_argument(
        "--low_complexity_filter",
        type=int,
        default=30, # default for fastp
        help="The threshold for low complexity filtering during quality control, a number between 0-100"
    )

    parser.add_argument(
        "--adapters",
        type=str,
        default="default_adapters.fa", 
        help="The path to your adapter file for adapter filtering"
    )

    parser.add_argument(
        "--snap_human_index_location",
        type=str,
        help="Location for the SNAP HG38 database index "
    )

    parser.add_argument(
        "--snap_edit_distance",
        type=int,
        default=27,
        help="The largest edit distance to keep for matches"
    )

    # below 4 flags are redundant assuming that a clinical computer should already have this built
    # for SNAP
    # flexible as e.g. it can build the entire NCBI NT, or Human HG38 genome, etd
    parser.add_argument(
        "--build_index",
        action="store_true"
        help="Builds an index from the genomes provided"
    )

    parser.add_argument(
        "--build_index_genomes",
        type=str,
        help="The genomes that you want to build an index for"
    )

    parser.add_argument(
        "--build_index_location",
        type=str,
        default="new_index"
        help="Builds an index from the genomes provided"
    )

    parser.add_argument(
        "--build_index_seed",
        type=int,
        default=27, #default seed for SNAP
        help="Seed to use for snap index"
    )

    # alternate to above 4:
    parser.add_argument(
        "--snap_index_location",
        type=str,
        help="Location for the SNAP NCBI NT database index"
    )

    # seed = portion of read to index (the larger the better but at a perforamnce cost)
    parser.add_argument(
        "--snap_index_seed",
        type=int,
        default=25, #default seed for SNAP
        help="Seed to use for snap index"
    )

    # The edit distancebetween two strings is the smallest number of single-character changes, insertions or deletions that are needed to change one into the other
    parser.add_argument(
        "--snap_index_edit_distance",
        type=int,
        default=25, #default seed for SNAP
        help="Seed to use for snap index"
    )

    # also, I do not believe that I will need index building commands for DIAMOND as I assume that they should already be there
    # diamond aligning
    parser.add_argument(
        "--diamond_index_location",
        type=str,
        help="Location for DIAMOND index"
    )

    # might be useful, might be not, more info in diamond man
    parser.add_argument(
        "--diamond_min-orf_len",
        type=int,
        default=30, # default for --min
        help="Location for DIAMOND index, Ignore translated sequences that do not contain an open reading frame of at least this length. Set to 1 to disable this."
    )

    parser.add_argument(
        "--diamond_e_value",
        type=int,
        default=-1, # default for --min
        help="E-value cutoff for diamond. Set to -1 to ignore."
    )

    # should we have options for frameshift, gapopen, gapextend penalties, and scoring matrix
    # also an option for tantan repeat masking, composition matrix

    # will we have access to a taxonmap/taxon names for make db


    parser_detect.set_defaults(func=run_pipeline)
    return parser.parse_args()


def main():
    
    # centrifuge is not installed locally
    path_command = subprocess.run('which centrifuge', 
                                    shell=True, 
                                    capture_output=True
                                )

    if (path_command.returncode != 0):
        raise Exception("centrifuge is not installed or initialised")

    args = parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
