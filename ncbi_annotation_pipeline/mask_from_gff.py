import argparse

import pandas as pd


def convert_to_sites(df: pd.DataFrame) -> set[int]:
    sites: set[int] = set()
    for start, end in zip(df["start"], df["end"]):
        # gff is 1-indexed inclusive
        sites.update(range(start - 1, end))
    return sites


def sites_to_file(sites: set[int], filename: str) -> None:
    with open(filename, "w", encoding="utf-8") as f:
        for site in sorted(sites):
            f.write(f"{site}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert GFF3 file to mask format for repeat regions and PE/PPE genes."
    )
    parser.add_argument(
        "gff_file",
        type=str,
        help="Input GFF3 file.",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Output file for the mask.",
    )
    args = parser.parse_args()

    gff = pd.read_csv(
        args.gff_file,
        sep="\t",
        header=None,
        comment="#",
        names=[
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )

    repeats = gff[gff["feature"] == "direct_repeat"]

    sites_to_file(
        convert_to_sites(repeats),
        args.output_file,
    )


if __name__ == "__main__":
    main()
