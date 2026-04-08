import argparse

import pandas as pd

from masquerade.expand_mask import enclose_mask, expand_mask


def report_status(df):
    """Report the status of the mask evaluation."""

    snp_locations = len(df[df["count"] > 0])
    high_frequency_snps = len(df[df["count"] > 5])
    very_high_frequency_snps = len(df[df["count"] > 40])
    ebr_locations = len(df[df["EBR"] < 1])
    low_ebr_locations = len(df[df["EBR"] < 0.9])
    very_low_ebr_locations = len(df[df["EBR"] < 0.1])

    print(f"Total positions: {len(df)}")
    print(
        f"SNP locations: {snp_locations} of which {high_frequency_snps} have high frequency (>5) and {very_high_frequency_snps} have very high frequency (>40)"
    )
    print(
        f"EBR locations: {ebr_locations} of which {low_ebr_locations} have EBR < 0.9 and {very_low_ebr_locations} have EBR < 0.1"
    )
    print("")

    return {
        "ont_snps": snp_locations,
        "ont_snps_high_frequency": high_frequency_snps,
        "ont_snps_very_high_frequency": very_high_frequency_snps,
        "ebr": ebr_locations,
        "ebr_low": low_ebr_locations,
        "ebr_very_low": very_low_ebr_locations,
    }



def main():
    parser = argparse.ArgumentParser(
        description="Evaluate mask sequence against SNP density."
    )
    parser.add_argument(
        "snp_density_file", type=str, help="Path to the SNP density csv."
    )
    parser.add_argument("ebr_file", type=str, help="Path to the EBR score csv.")
    # parser.add_argument("illumina_false_snps_file", type=str, help="Path to the Illumina false SNPs mask according to marin.")
    parser.add_argument(
        "-o",
        "--output_root",
        default="evaluation",
        type=str,
        help="Path to save the output files.",
    )
    parser.add_argument("mask", nargs="+", help="List of mask files to process")

    args = parser.parse_args()
    output_root = args.output_root

    # Read SNP density data
    snp_density = pd.read_csv(args.snp_density_file)
    ebr_df = pd.read_csv(args.ebr_file)
    # illumina_false_snps = pd.read_csv(args.illumina_false_snps_file, header=None, names=["Pos"])
    # illumina_false_snps["illumina_false_snp"] = True

    df = snp_density.merge(ebr_df, on="Pos", how="outer").fillna({"EBR": 1, "count": 0})
    print("starting evaluation")

    current_mask = set()
    rows = []

    first_row = {
        "mask": "initial",
        "mask_size": len(current_mask),
        "additional_masking": 0,
        "total_mask_size": len(current_mask),
        "type": "current_status",
        **report_status(df),
    }
    rows.append(first_row)

    for mask_file in args.mask:
        if mask_file.startswith("expand"):
            mask_name = mask_file
            amount = int(mask_name.replace("expand", ""))
            mask = expand_mask(current_mask, amount)
        elif mask_file.startswith("enclose"):
            mask_name = mask_file
            amount = int(mask_name.replace("enclose", ""))
            mask = enclose_mask(current_mask, amount)
        else:
            mask_name = mask_file.split("/")[-1].split(".")[0]
            mask = set(pd.read_csv(mask_file, header=None, names=["Pos"])["Pos"])
        new_sites = mask - current_mask
        current_mask.update(mask)
        print(
            f"### Processing mask: {mask_name} with {len(mask)} positions of which {len(new_sites)} are new."
        )

        print(f"Current mask size: {len(current_mask)} positions.")
        print("This leaves the following positions unmasked:")
        base_data = {
            "mask": mask_name,
            "mask_size": len(mask),
            "additional_masking": len(new_sites),
            "total_mask_size": len(current_mask),
        }

        data_change = {
            **base_data,
            "type": "extra_masking",
            **report_status(df[df["Pos"].isin(new_sites)]),
        }

        data_current = {
            **base_data,
            "type": "current_status",
            **report_status(df[~df["Pos"].isin(current_mask)]),
        }

        rows.append(data_change)
        rows.append(data_current)

    # Save the results to a CSV file
    output_df = pd.DataFrame(rows)

    output_df[output_df["type"] == "current_status"].to_csv(
        f"{output_root}.csv", index=False
    )
    output_df[output_df["type"] == "extra_masking"].to_csv(
        f"{output_root}_extra_masking.csv", index=False
    )


if __name__ == "__main__":
    main()
