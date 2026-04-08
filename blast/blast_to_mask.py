import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Convert BLAST output to mask format.")
    parser.add_argument("blast_file", type=str, help="Input BLAST output file.")
    parser.add_argument("mask_file", type=str, help="Output mask file.")
    parser.add_argument(
        "--min_length", type=int, default=0, help="Minimum length of region."
    )
    parser.add_argument(
        "--min_identity", type=float, default=0, help="Minimum identity percentage."
    )
    parser.add_argument(
        "--max_evalue", type=float, default=1, help="Minimum e-value."
    )

    args = parser.parse_args()

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
    df = df.query('qstart != sstart or qend != send')
    
    df = df.query(
        "length >= @args.min_length and pident >= @args.min_identity and evalue <= @args.max_evalue"
    )
    
    sites = set()
    for _index, row in df.iterrows():
        for site in range(row["qstart"], row["qend"] + 1):
            sites.add(site)
        for site in range(row["sstart"], row["send"] + 1):
            sites.add(site)
    print(df)
    
    ordered_sites = sorted(list(sites))
    with open(args.mask_file, "w", encoding="utf-8") as f:
        for site in ordered_sites:
            # output as 0-indexed to match the mask format
            f.write(f"{site - 1}\n")
    print(f"Masked {len(ordered_sites)} sites.")


if __name__ == "__main__":
    main()
