import argparse
import sys 
from main_pipeline import detection

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
