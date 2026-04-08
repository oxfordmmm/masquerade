import shutil
import subprocess
import os
from pathlib import Path
import argparse
import pandas as pd

from masquerade.expand_mask import enclose_mask

def run_self_blast(fasta_path, threads=4):
    fasta_path = Path(fasta_path)
    db_dir = Path("tmp_blast_db")
    db_prefix = db_dir / "db"
    outfile = "tmp_blast_results.tsv"

    # Ensure output directory exists
    db_dir.mkdir(parents=True, exist_ok=True)

    # 1. Make BLAST database
    makeblastdb_cmd = [
        "makeblastdb",
        "-in",
        str(fasta_path),
        "-dbtype",
        "nucl",
        "-out",
        str(db_prefix),
        "-title",
        "DB",
    ]
    subprocess.run(makeblastdb_cmd, check=True)

    # 2. Run blastn
    blastn_cmd = [
        "blastn",
        "-task",
        "megablast",
        "-query",
        str(fasta_path),
        "-db",
        str(db_prefix),
        "-dust",
        "no",
        "-outfmt",
        "6 qseqid qstart qend sseqid sstart send length mismatch pident evalue",
        "-evalue",
        "0.0001",
        "-num_threads",
        str(threads),
        "-word_size",
        "15",
        "-gapopen",
        "5",
        "-gapextend",
        "2",
        "-max_target_seqs",
        "10000000",
        "-out",
        outfile,
    ]
    subprocess.run(blastn_cmd, check=True)

    # 3. Load into pandas
    cols = [
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
    ]

    df = pd.read_csv(outfile, sep="\t", names=cols)

    # remove temporary files
    os.remove(outfile)
    shutil.rmtree(db_dir)
    

    return df


def blast_hits_to_sites(df, min_length=100, min_identity=90.0) -> set[int]:
    # Remove self-hits (diagonal)
    df = df.query("qstart != sstart or qend != send")

    # Filter based on criteria
    df = df[(df["length"] >= min_length) & (df["pident"] >= min_identity)]

    sites = set()
    for _index, row in df.iterrows():
        for site in range(row["qstart"], row["qend"] + 1):
            sites.add(site)
        for site in range(row["sstart"], row["send"] + 1):
            sites.add(site)

    return sites


def run_dustmasker(fasta_path, level: int) -> set[int]:
    fasta_path = Path(fasta_path)
    outfile = "tmp_dust.tsv"

    # 1. Run dustmasker
    cmd = [
        "dustmasker",
        "-in",
        str(fasta_path),
        "-infmt",
        "fasta",
        "-parse_seqids",
        "-outfmt",
        "interval",
        "-out",
        outfile,
        "-level",
        str(level),
    ]
    subprocess.run(cmd, check=True)

    # 2. Parse interval output
    df = pd.read_csv(
        outfile,
        sep=" - ",
        header=None,
        names=[
            "start",
            "end",
        ],
        skiprows=1,
        engine="python",
    )

    sites = set()
    for start, end in zip(df["start"], df["end"]):
        # is inclusive
        for site in range(start, end + 1):
            sites.add(site)

    os.remove(outfile)
    return sites


def run_mummer_tandems(fasta_path, min_match=31, min_extent=75) -> set[int]:
    fasta_path = Path(fasta_path)
    outfile = "tmp_mummer_tandems.tsv"

    # 1. Run exact-tandems
    cmd = f"exact-tandems {fasta_path} {min_match} > {outfile}"
    subprocess.run(cmd, shell=True, check=True)

    # 2. Read output into DataFrame
    tandem_df = pd.read_csv(outfile, sep=r"\s+")

    sites = set()
    for start, extent in zip(tandem_df["Start"], tandem_df["Extent"]):
        if extent < min_extent:
            continue
        for site in range(start, start + extent):
            sites.add(site)

    os.remove(outfile)
    return sites


def main():
    parser = argparse.ArgumentParser(description="Auto Masking Script")
    parser.add_argument(
        "fasta", type=str, help="Path to fasta file to mask"
    )
    parser.add_argument(
        "output", type=str, help="Output file name"
    )
    parser.add_argument(
        "--min_length",
        type=int,
        default=100,
        help="Minimum length of BLAST hit to consider",
    )
    parser.add_argument(
        "--min_identity",
        type=float,
        default=90.0,
        help="Minimum identity percentage of BLAST hit to consider",
    )
    parser.add_argument(
        "--threads", type=int, default=4, help="Number of threads to use for BLAST"
    )
    parser.add_argument(
        "--dust_level",
        type=int,
        default=20,
        help="Dust level (0-20) for dustmasker (higher is more aggressive)",
    )
    parser.add_argument(
        "--min_tandem_extent",
        type=int,
        default=75,
        help="Minimum extent of tandem repeat required to mask it",
    )
    parser.add_argument(
        "--min_tandem_match",
        type=int,
        default=31,
        help="Minimum match length for exact-tandems",
    )
    parser.add_argument(
        "--one_indexed",
        action="store_true",
        help="Output 1-indexed masked sites instead of 0-indexed",
    )
    parser.add_argument(
        "--enclose",
        type=int,
        default=100,
        help="Enclose gaps in the mask of size <= this value",
    )
    args = parser.parse_args()

    print("Running self-BLAST...")
    blast_df = run_self_blast(args.fasta, threads=args.threads)
    print("Processing BLAST hits...")
    blast_mask = blast_hits_to_sites(
        blast_df, min_length=args.min_length, min_identity=args.min_identity
    )

    print("Running dustmasker...")
    dust_mask = run_dustmasker(args.fasta, level=args.dust_level)

    print("Running mummer exact-tandems...")
    tandem_mask = run_mummer_tandems(
        args.fasta, min_match=args.min_tandem_match, min_extent=args.min_tandem_extent
    )

    print("Combining masks and enclosing gaps...")
    all_masked_sites = blast_mask.union(dust_mask).union(tandem_mask)
    if args.enclose > 0:
        all_masked_sites = enclose_mask(all_masked_sites, args.enclose)

    ordered_sites = sorted(list(all_masked_sites))
    with open(args.output, "w", encoding="utf-8") as f:
        if not args.one_indexed:
            ordered_sites = [site - 1 for site in ordered_sites]

        for site in ordered_sites:
            f.write(f"{site}\n")
    print(f"Masked {len(ordered_sites)} sites.")
