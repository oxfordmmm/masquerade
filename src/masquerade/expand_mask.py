import argparse
import sys

import pandas as pd


def expand_mask(mask: set[int], amount: int) -> set[int]:
    """Expand the mask by a given amount on both sides."""

    expanded_mask = mask.copy()
    for pos in mask:
        for direction in [-1, 1]:
            for i in range(1, amount + 1):
                new_pos = pos + direction * i
                if new_pos not in expanded_mask:
                    expanded_mask.add(new_pos)
                else:
                    break  # Stop expanding in this direction if we hit an existing position

    return expanded_mask


def enclose_mask(mask: set[int], amount: int) -> set[int]:
    """Any gaps in the mask of size <=amount will be filled in."""
    expanded_mask = mask.copy()
    sorted_mask = sorted(mask)

    for i in range(len(sorted_mask) - 1):
        start = sorted_mask[i]
        end = sorted_mask[i + 1]

        if 0 < end - start - 1 <= amount:
            expanded_mask.update(range(start + 1, end))

    return expanded_mask


def main():
    parser = argparse.ArgumentParser(
        description="Compare two TSV files and report differences."
    )
    parser.add_argument("input", type=str, help="Path to the input file.")
    parser.add_argument(
        "-a", "--amount", type=int, default=0, help="Amount to expand mask."
    )
    parser.add_argument(
        "-e",
        "--enclose",
        type=int,
        default=0,
        help="Enclose gaps in the mask if set to a positive integer.",
    )
    parser.add_argument("output", type=str, help="Path to the output file.")
    args = parser.parse_args()

    if args.input == "-":
        df = pd.read_csv(sys.stdin, sep="\t", header=None, names=["site"])
    else:
        df = pd.read_csv(args.input, sep="\t", header=None, names=["site"])

    # Load the TSV files
    mask = set(df["site"])
    expanded_mask = mask.copy()
    if args.amount > 0:
        expanded_mask = expand_mask(mask, args.amount)
    if args.enclose > 0:
        expanded_mask = enclose_mask(expanded_mask, args.enclose)
    expanded_df = pd.DataFrame(sorted(expanded_mask), columns=["site"])
    expanded_df.to_csv(args.output, sep="\t", index=False, header=False)

    print("Expanded mask sizes:", len(mask), "->", len(expanded_mask))


if __name__ == "__main__":
    main()
