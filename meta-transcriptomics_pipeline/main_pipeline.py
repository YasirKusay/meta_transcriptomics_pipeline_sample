import argparse

def run_pipeline(args: argparse.Namespace):
    results = detect_contamination(args)
    if (args.graphical_output):
        draw(results, args.output_file)
    else:
        table(results)
