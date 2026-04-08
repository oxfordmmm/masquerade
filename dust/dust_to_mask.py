import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Convert dust output to mask format.")
    parser.add_argument("dust_file", type=str, help="Input dust output file.")
    parser.add_argument("mask_file", type=str, help="Output mask file.")

    args = parser.parse_args()

    col_names = [
        "start",
        "end",
    ]
    df = pd.read_csv(args.dust_file, sep=" - ", header=None, names=col_names, skiprows=1)
    print(df)
    
    sites = set()
    for start, end in zip(df["start"], df["end"]):
        # is inclusive
        for site in range(start, end + 1):
            sites.add(site)
    
    ordered_sites = sorted(list(sites))
    with open(args.mask_file, "w", encoding="utf-8") as f:
        for site in ordered_sites:
            # output as 0-indexed to match the mask format
            f.write(f"{site - 1}\n")
    print(f"Masked {len(ordered_sites)} sites.")


if __name__ == "__main__":
    main()
