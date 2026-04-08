import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Find all rows in blastn output containing given site."
    )
    parser.add_argument("blast_file", type=str, help="Input BLAST output file.")
    parser.add_argument("site", type=int, help="Site to find.")

    args = parser.parse_args()

    site = int(args.site)

    col_names = [
        "qseqid",
        "qstart",
        "qend",
        "sseqid",
        "sstart",
        "send",
        "length",
        "mismatch",
        "pident",
        "evalue",
    ]
    df = pd.read_csv(args.blast_file, sep="\t", header=None, names=col_names)

    # remove diagonal hits
    df = df.query("qstart != sstart or qend != send")

    query_matches = df.query("qstart <= @site and qend >= @site")
    subject_matches = df.query("sstart <= @site and send >= @site")
    if query_matches.empty and subject_matches.empty:
        print(f"No matches found for site {site}.")
        return

    combined_matches = pd.concat([query_matches, subject_matches]).drop_duplicates()
    print(f"Matches found for site {site}:")
    print(combined_matches)


if __name__ == "__main__":
    main()
