# Takes a fasta density summary, and some masks and marks which masks the sites are in

import argparse

import pandas as pd


def find_distance_to_mask(pos: int, mask: set[int]) -> int:
    if pos in mask:
        return 0

    return min(abs(pos - m) for m in mask)


def find_directional_distances_to_mask(pos: int, mask: set[int]) -> tuple[int, int]:
    """
    Find the distance to the nearest mask position in both directions.
    Returns a tuple of (distance to left, distance to right).
    If there is no mask position in one direction, returns None for that direction.
    """
    if pos in mask:
        return 0, 0

    left_distance = min(((pos - m) for m in mask if m < pos), default=999999)
    right_distance = min(((m - pos) for m in mask if m > pos), default=999999)

    return (left_distance, right_distance)


def main():
    parser = argparse.ArgumentParser(
        description="Look up locations in a fasta density summary and mark which masks they are in."
    )
    parser.add_argument(
        "fasta_density_summary",
        type=str,
        help="Path to the fasta density summary file.",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to the output file. This will be a csv file with the same format as the input file.",
    )
    parser.add_argument(
        "--find_distance",
        action="store_true",
        help="If set, will find the distance to the nearest mask for each position.",
    )
    parser.add_argument(
        "--find_dir_distance",
        action="store_true",
        help="If set, will find the distance to the nearest mask in both directions for each position.",
    )
    parser.add_argument(
        "masks",
        nargs="+",
        type=str,
        help="Path to the masks file. Expected to be 0-indexed.",
    )

    args = parser.parse_args()

    density = pd.read_csv(args.fasta_density_summary, usecols=["Pos", "count"])
    mask_cols = []

    for mask in args.masks:
        mask_name = mask.split("/")[-1].split(".")[0]
        mask_cols.append(mask_name)
        mask_df = pd.read_csv(mask, header=None, names=["pos"])
        mask_set = set(mask_df["pos"])

        if args.find_distance:
            density[mask_name] = density["Pos"].apply(
                lambda pos: find_distance_to_mask(pos, mask_set)
            )
        elif args.find_dir_distance:
            density[mask_name + "_left"], density[mask_name + "_right"] = zip(
                *density["Pos"].apply(lambda pos: find_directional_distances_to_mask(pos, mask_set))
            )
        else:
            # positions from compare_fasta are already 0-indexed like the masks
            density[mask_name] = density["Pos"].isin(mask_set)

    if args.find_distance:
        density["closest_mask"] = density[mask_cols].min(axis=1)

    print(density)
    density.to_csv(args.output, index=False)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
