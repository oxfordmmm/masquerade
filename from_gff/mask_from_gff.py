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
        default="../h37rv.gff3",
        help="Input GFF3 file.",
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

    gff["gene"] = gff["attribute"].str.extract(r"gene=([^;]+)")
    gff["gene_synonym"] = gff["attribute"].str.extract(r"gene_synonym=([^;]+)")
    print(gff)

    repeat_regions = gff[gff["feature"] == "repeat_region"][["start", "end"]]
    repeat_regions["type"] = "repeat_region"
    repeat_regions["gene"] = ""
    sites_to_file(
        convert_to_sites(repeat_regions),
        "repeat_regions.mask",
    )

    ppe_pe = gff[
        (gff["feature"] == "gene")
        & (
            (gff["gene"].str.contains("PPE|PE"))
            | (gff["gene_synonym"].str.contains("PPE|PE"))
        )
    ]
    ppe_pe = ppe_pe[["start", "end", "gene"]]
    ppe_pe["type"] = "PE/PPE"
    sites_to_file(
        convert_to_sites(ppe_pe),
        "ppe_pe.mask",
    )

    regions = pd.concat([repeat_regions, ppe_pe], ignore_index=True)
    regions = regions[["start", "end", "type", "gene"]]
    regions = regions.sort_values(by=["start", "end"])
    regions.to_csv("regions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
