import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Convert mummer output to mask format.")
    parser.add_argument("--input", type=str, help="Input mummer output file.")
    parser.add_argument("--tandem_file", type=str, help="Input tandem file.")
    parser.add_argument("--min_extent", type=int, default=75, help="Minimum extent of tandem repeat.")
    parser.add_argument("mask_file", type=str, help="Output mask file.")
    parser.add_argument(
        "--min_length", type=int, default=75, help="Minimum length of region."
    )
    parser.add_argument(
        "--min_identity", type=float, default=90, help="Minimum identity percentage."
    )

    args = parser.parse_args()

    col_names = [
        "qstart",
        "qend",
        "sstart",
        "send",
        "aln_length1",
        "aln_length2",
        "pident",
        "length1",
        "length2",
        "cov1",
        "cov2",
        "tag1",
        "tag2",
    ]
    
    sites = set()
    if args.input:
        
        df = pd.read_csv(args.input, sep="\t", header=None, skiprows=4, names=col_names)
        print(df)
        # remove diagonal hits
        df = df.query("qstart != sstart or qend != send")
        df = df.query("aln_length1 >= @args.min_length and pident >= @args.min_identity")
        
        for _index, row in df.iterrows():
            for site in range(row["qstart"], row["qend"] + 1):
                sites.add(site)
            for site in range(row["sstart"], row["send"] + 1):
                sites.add(site)
        print(df)
    
    
    if args.tandem_file:
        tandem_df = pd.read_csv(args.tandem_file, sep=r'\s+')
        print(tandem_df)
        for start, extent in zip(tandem_df["Start"], tandem_df["Extent"]):
            if extent < args.min_extent:
                continue
            for site in range(start, start + extent):
                sites.add(site)

    ordered_sites = sorted(list(sites))
    with open(args.mask_file, "w", encoding="utf-8") as f:
        for site in ordered_sites:
            # output as 0-indexed to match the mask format
            f.write(f"{site - 1}\n")
    print(f"Masked {len(ordered_sites)} sites.")


if __name__ == "__main__":
    main()
