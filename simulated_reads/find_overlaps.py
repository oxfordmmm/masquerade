import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Find overlaps.")
    parser.add_argument("blast_results", type=str, help="Path to the blast results file.")
    parser.add_argument(
        "-m",
        "--mask-output",
        type=str,
        default="overlaps.mask",
        help="Path to the output file for overlaps mask.",
    )
    parser.add_argument(
        "-b",
        "--blast-output",
        type=str,
        default="overlaps.tsv",
        help="Path to the output file for blast results.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=50,
        help="Minimum length of the overlap to consider.",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=90,
        help="Minimum identity of the overlap to consider.",
    )
    args = parser.parse_args()

    blast_df = pd.read_csv(
        args.blast_results,
        sep="\t",
        header=None,
        names=[
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
        ],
    )

    # remove diagonal mapping
    blast_df = blast_df[blast_df["qseqid"] != blast_df["sstart"]]

    # drop duplicate targets
    blast_df = blast_df.drop_duplicates(subset=["sstart", "send"], keep="first")

    # filter by length and identity
    blast_df = blast_df[blast_df["length"] >= args.min_length]
    blast_df = blast_df[blast_df["pident"] >= args.min_identity]

    sites: set[int] = set()
    for qseqid, qstart, qend, sstart, send in zip(
        blast_df["qseqid"],
        blast_df["qstart"],
        blast_df["qend"],
        blast_df["sstart"],
        blast_df["send"],
    ):
        sites.update(
            range(
                qseqid + qstart - 1,
                qseqid + qend,
            )
        )

        sites.update(
            range(
                sstart,
                send + 1,
            )
        )

    blast_df.to_csv(
        args.blast_output, index=False, sep="\t" if args.blast_output.endswith(".tsv") else ","
    )

    with open(args.mask_output, "w", encoding="utf-8") as mask_file:
        for i in sorted(sites):
            # convert to 0-based index
            mask_file.write(f"{i-1}\n")
    
    print(f"Mask file written to {args.mask_output} contains {len(sites)} sites.")


if __name__ == "__main__":
    main()
