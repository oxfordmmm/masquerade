import argparse


def read_matches(alignments: str) -> list[tuple[int, int]]:
    """
    Read matches from a string and return a list of tuples (seq, pos).
    """
    if alignments == "":
        return []
    matches = [
        (int(match.split(",")[0]), int(match.split(",")[1]))
        for match in alignments.split("|")
    ]

    return matches


def main():
    parser = argparse.ArgumentParser(description="Create a mask from genmap csv.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input CSV file")
    parser.add_argument(
        "-k", "--kmer_size", type=int, default=50, help="K-mer size for expansion"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Output mask file"
    )

    args = parser.parse_args()
    kmer_size = args.kmer_size
    mask = set()

    with open(args.input, "r") as infile:
        for line in infile:
            if line.startswith('"'):
                # this is first line
                continue
            query, pos_strand, neg_strand = line.strip().split(";")
            seq, pos = map(int, query.split(","))

            matches = read_matches(pos_strand) + read_matches(neg_strand)

            if len(matches) > 1:
                for i in range(kmer_size):
                    mask.add(pos + i)

    with open(args.output, "w") as outfile:
        for pos in sorted(mask):
            outfile.write(f"{pos}\n")


if __name__ == "__main__":
    main()
