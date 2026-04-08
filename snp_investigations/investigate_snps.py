import argparse
import os
import subprocess

import pandas as pd
from Bio import SeqIO

# masks, and requested regions assumed to be 0-indexed


def parse_range(range_str):
    if ":" not in range_str:
        # Then it is a single position so do range around point
        middle = int(range_str)
        start = middle - 200
        end = middle + 200
    else:
        start, end = map(int, range_str.split(":"))
    return start, end


def read_snps_csv(snps_csv, start, end) -> pd.DataFrame:
    df = pd.read_csv(snps_csv)
    # Filter SNPs within the specified range
    filtered_df = df[(df["Pos"] >= start) & (df["Pos"] <= end)]
    return filtered_df


def read_masks(masks: list[str], start, end) -> pd.DataFrame:
    dfs = []
    for mask_path in masks:
        mask_name = mask_path.split("/")[-1].split(".")[0]
        df = pd.read_csv(mask_path, header=None, names=["pos"])
        df = df[(df["pos"] >= start) & (df["pos"] <= end)]
        df["mask_name"] = mask_name
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)

    # pivot wide
    result = df.pivot_table(
        index="pos", columns="mask_name", aggfunc=lambda x: True, fill_value=False
    )
    result.reset_index(inplace=True)
    result.rename(columns={"pos": "Pos"}, inplace=True)

    return result


def get_subsequence(fasta_file, start, end):
    record = SeqIO.read(fasta_file, "fasta")
    return str(record.seq[start : end + 1])


def get_msa_from_dir(fasta_dir, samples, start, end) -> str:
    msa_str = ""
    for sample in samples:
        fasta_path = os.path.join(fasta_dir, f"{sample}.final.fasta")
        if os.path.exists(fasta_path):
            seq = get_subsequence(fasta_path, start, end)
            msa_str += f"{seq}\n"
        else:
            msa_str += "NOT_FOUND\n"
    return msa_str


def run_blast(query_fasta, blast_db, output_file, start):
    command = (
        f"blastn -task megablast -query {query_fasta} -db {blast_db} "
        f"-dust no -outfmt '6 qseqid qstart qend sseqid sstart send length mismatch pident evalue' "
        f"-evalue 0.0001 -num_threads 1 "
        f"-word_size 17 -gapopen 5 -gapextend 2 "
        f"-out {output_file}"
    )
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

    blast_df = pd.read_csv(
        output_file,
        sep="\t",
        header=None,
        names=[
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
        ],
    )
    blast_df["qstart"] += start
    blast_df["qend"] += start
    # filter out diagonal hits
    blast_df = blast_df.query('qstart != sstart or qend != send')
    
    # Save the filtered results
    blast_df.to_csv(output_file, sep="\t", index=False)
    print(f"BLAST results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Extract SNP region sequences from samples."
    )
    parser.add_argument(
        "--range", required=True, help="Range in format start:end (e.g., 100:500)"
    )
    parser.add_argument("--ref_fasta", required=True, help="Reference fasta file")
    parser.add_argument(
        "--snps_csv", required=True, help="SNPs CSV file (columns: sample,position,...)"
    )
    parser.add_argument(
        "--ont_fasta_dir", required=True, help="Directory of ont sample fasta files"
    )
    parser.add_argument(
        "--illumina_fasta_dir",
        required=True,
        help="Directory of illumina sample fasta files",
    )
    parser.add_argument(
        "--masks",
        help="Path to masks files",
        nargs="+",
        required=False,
    )
    parser.add_argument(
        "--blast_db",
        help="Path to BLAST database for SNPs",
        required=False,
    )
    parser.add_argument("--output", required=True, help="Output root for files")
    args = parser.parse_args()
    output_root = args.output

    start, end = parse_range(args.range)

    snps_df = read_snps_csv(args.snps_csv, start, end)
    snps_df.to_csv(output_root + "snps.csv", index=False)
    samples = snps_df["Sample"].unique().tolist()

    # Get reference sequence
    ref_seq = get_subsequence(args.ref_fasta, start, end)
    with open(output_root + "ref.fasta", "w", encoding="utf-8") as out:
        out.write(f">ref\n{ref_seq}\n")

    # Get mask info about the region
    if args.masks:
        masks_df = read_masks(args.masks, start, end)
        snp_count = snps_df.groupby("Pos").size().reset_index(name="SNP_count")
        if masks_df.empty:
            print("No masks found in the specified range.")
            snp_count.to_csv(output_root + "masking.csv", index=False)
        else:
            masks_df = pd.merge(snp_count, masks_df, on="Pos", how="outer").copy()
            # fill snp count NaNs with 0
            masks_df["SNP_count"] = masks_df["SNP_count"].fillna(0).astype(int)

            masks_df.fillna(False, inplace=True)
            masks_df["Pos"] = masks_df["Pos"].astype(int)
            masks_df = masks_df.sort_values(by="Pos").reset_index(drop=True)
            masks_df.to_csv(output_root + "masking.csv", index=False)

    # Run BLAST if a database is provided
    if args.blast_db:
        blast_output_file = output_root + "blast_results.tsv"
        run_blast(output_root + "ref.fasta", args.blast_db, blast_output_file, start)
        print(f"BLAST results saved to {blast_output_file}")

    # Want to mark sites which have snps
    snps = snps_df["Pos"].tolist()
    marker_str = [" "] * (end - start + 1)
    for pos in snps:
        if start <= pos <= end:
            marker_str[pos - start] = "X"

    with open(output_root + "ont.csv", "w", encoding="utf-8") as out:
        out.write("sample,sequence\n")
        out.write(f"{ref_seq},ref\n")
        out.write(f"{''.join(marker_str)},marker\n")

        ont_msa = get_msa_from_dir(args.ont_fasta_dir, samples, start, end)
        out.write(ont_msa)

    with open(output_root + "illumina.csv", "w", encoding="utf-8") as out:
        out.write("sample,sequence\n")
        out.write(f"{ref_seq},ref\n")
        out.write(f"{''.join(marker_str)},marker\n")
        illumina_msa = get_msa_from_dir(args.illumina_fasta_dir, samples, start, end)
        out.write(illumina_msa)


if __name__ == "__main__":
    main()
