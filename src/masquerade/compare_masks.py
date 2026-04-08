import argparse

import pandas as pd

from masquerade.expand_mask import expand_mask


def main():
    parser = argparse.ArgumentParser(
        description="Compare two TSV files and report differences."
    )
    parser.add_argument("file1", type=str, help="Path to the first TSV file.")
    parser.add_argument("file2", type=str, help="Path to the second TSV file.")
    parser.add_argument("--snps", type=str, help="Path to the SNP density file.")
    args = parser.parse_args()

    # Load the TSV files
    df1 = pd.read_csv(args.file1, sep="\t", header=None, names=["Pos"])
    df2 = pd.read_csv(args.file2, sep="\t", header=None, names=["Pos"])

    mask1 = set(df1["Pos"])
    mask2 = set(df2["Pos"])

    total_sites = len(mask1.union(mask2))

    expanded_mask1 = expand_mask(mask1, 50)
    expanded_mask2 = expand_mask(mask2, 50)
    expanded_df1 = pd.DataFrame(sorted(expanded_mask1), columns=["Pos"])
    expanded_df2 = pd.DataFrame(sorted(expanded_mask2), columns=["Pos"])

    df1["mask1"] = True
    df2["mask2"] = True
    expanded_df1["expanded_mask1"] = True
    expanded_df2["expanded_mask2"] = True

    # Merge the dataframes to find differences
    merged_df = (
        pd.merge(df1, df2, on="Pos", how="outer")
        .merge(expanded_df1, on="Pos", how="outer")
        .merge(expanded_df2, on="Pos", how="outer")
        .fillna(False)
    )

    merged_df.sort_values(by="Pos", inplace=True)

    # only want sites where masks differ
    merged_df = merged_df[merged_df["mask1"] != merged_df["mask2"]]
    merged_df["is_close"] = merged_df["expanded_mask1"] & merged_df["expanded_mask2"]

    unique_to_1 = len(merged_df[merged_df["mask1"]])
    unique_to_2 = len(merged_df[merged_df["mask2"]])

    unique_to_1_but_close_to_2 = len(
        merged_df[merged_df["mask1"] & merged_df["is_close"]]
    )
    unique_to_2_but_close_to_1 = len(
        merged_df[merged_df["mask2"] & merged_df["is_close"]]
    )

    print(f"Total unique sites in both files: {total_sites}")
    print(f"Sites in both masks: {len(mask1.intersection(mask2))}")
    print(f"Unique to file 1: {unique_to_1}")
    print(f"Unique to file 2: {unique_to_2}")
    print(f"Unique to file 1 but close to file 2: {unique_to_1_but_close_to_2}")
    print(f"Unique to file 2 but close to file 1: {unique_to_2_but_close_to_1}")

    if args.snps:
        snp_df = pd.read_csv(args.snps, usecols=["Pos", "count"]).rename(
            columns={"count": "snp_count"}
        )
        merged_df = pd.merge(merged_df, snp_df, on="Pos", how="left")
        merged_df.fillna({"snp_count": 0}, inplace=True)

    merged_df.to_csv("comparison_output.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
