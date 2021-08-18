import argparse
import sys 
from detect_contamination import detection

def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect contamination in an assembly"
    )

    parser.add_argument(
        "--fragment_length",
        type=int,
        default=150,
        help="the length of the fragment sequence"
    )

    return parser.parse_args()


def main():
    args = parse_args()
    detection(args)

if __name__ == "__main__":
    main()
