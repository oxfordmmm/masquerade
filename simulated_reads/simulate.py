import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Simulate reads from a reference genome."
    )
    parser.add_argument("reference", type=str, help="Path to the reference genome file.")
    parser.add_argument(
        "output", type=str, help="Path to the output file for simulated reads."
    )
    parser.add_argument(
        "--copies", type=int, default=1, help="Number of copies of each read"
    )
    parser.add_argument(
        "--read-length", type=int, default=150, help="Length of each read."
    )
    parser.add_argument("--spacing", type=int, default=1, help="Spacing between reads.")

    args = parser.parse_args()
    copies = int(args.copies)

    # read the reference genome
    with open(args.reference, "r", encoding="utf-8") as ref_file:
        lines = ref_file.readlines()
        fasta = "".join(line.strip() for line in lines if not line.startswith(">"))

    fasta_length = len(fasta)
    print(f"Fasta is {fasta_length}bp long")

    fasta = fasta + fasta[: args.read_length]

    with open(args.output, "w", encoding="utf-8") as out_file:
        i = 0
        while i < fasta_length:
            read = fasta[i : i + args.read_length]

            if copies > 1:
                for j in range(copies):
                    out_file.write(f">{i+1}:{j}\n")
                    out_file.write(f"{read}\n")
            else:
                out_file.write(f">{i+1}\n")
                out_file.write(f"{read}\n")

            i += args.spacing


if __name__ == "__main__":
    main()
