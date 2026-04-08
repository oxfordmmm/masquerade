import argparse
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Combine and deduplicate sites from multiple files."
    )
    parser.add_argument("files", nargs="+", help="List of input files to process")
    parser.add_argument(
        "-o", "--output", type=str, default="-", help="Output file name"
    )
    args = parser.parse_args()

    combined_sites = set()

    for file in args.files:
        try:
            data = pd.read_csv(file, header=None, names=["site"])
            combined_sites.update(data["site"].unique())
        except Exception as e:
            print(f"Error reading file {file}: {e}")

    ordered_sites = sorted(list(combined_sites))
    
    # If output is "-", print to stdout
    if args.output == "-":
        output_stream = sys.stdout
    else:
        output_stream = open(args.output, "w", encoding="utf-8")

    # Write the ordered sites to the output
    for site in ordered_sites:
        output_stream.write(f"{site}\n")

    if args.output != "-":
        output_stream.close()

    # write debug stats to stderr
    print(f"Combined {len(combined_sites)} unique sites from {len(args.files)} files.", file=sys.stderr)


if __name__ == "__main__":
    main()
