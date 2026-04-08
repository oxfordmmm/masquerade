import pandas as pd


def write_mask(mask: set[int], filename: str):
    with open(filename, "w") as f:
        for pos in sorted(mask):
            f.write(f"{pos}\n")


def read_bam(df: pd.DataFrame) -> set[int]:
    mask: set[int] = set()
    for start, end in zip(df["chromStart"], df["chromEnd"]):
        # bed files 0-indexed and end is exclusive
        mask.update(range(start, end))
    return mask


def read_ambiguous_regions():
    df = pd.read_csv(
        "supplement/AF5_CommonAmbigousRegions_H37Rv.WithHeader.bed.tsv", sep="\t"
    )

    write_mask(read_bam(df), "marin_ambiguous_regions.mask")


def refinded_low_confidence_regions():
    df = pd.read_csv(
        "supplement/AF13_RLC_Regions.H37Rv.bed",
        sep="\t",
        header=None,
        names=["chrom", "chromStart", "chromEnd"],
    )
    mask = read_bam(df)
    write_mask(mask, "marin_RLC.mask")


def snp_hotspots():
    df = pd.read_csv(
        "supplement/AF21_210202_Mtb_H37rv.Top30SourcesOfFPs.FiltMQ30.bed",
        sep="\t",
        header=None,
        names=[
            "chrom",
            "chromStart",
            "chromEnd",
            "strand",
            "H37rv_GeneID",
            "Symbol",
            "ExcludedGroup_Category",
            "PEandPPE_Subfamily",
            "Functional_Category",
            "FP_Count",
            "Length",
        ],
    )
    mask = read_bam(df)
    write_mask(mask, "marin_snp_hotspots.mask")

def empiral_base_level_recall():
    df = pd.read_csv(
        "supplement/AF18_H37Rv_EBR_36CI.bedgraph",
        sep="\t",
        header=None,
        names=["chrom", "chromStart", "chromEnd", "EBR"],
    )
    non_ambiguous_df = df[df["EBR"] != -1]  # remove the -1 values which are ambiguous
    
    write_mask(read_bam(non_ambiguous_df[non_ambiguous_df["EBR"] < 0.9]), "marin_EBR_sub_90.mask")
    write_mask(read_bam(non_ambiguous_df[non_ambiguous_df["EBR"] < 0.95]), "marin_EBR_sub_95.mask")
    write_mask(read_bam(non_ambiguous_df[non_ambiguous_df["EBR"] < 1]), "marin_EBR_sub_100.mask")
    
    df = df[df["EBR"] < 1]
    sites = []
    for start, end, ebr in zip(df["chromStart"], df["chromEnd"], df["EBR"]):
        # start in inclusive, end is exclusive
        for pos in range(start, end):
            sites.append((pos, ebr))

    mask_df = pd.DataFrame(sites, columns=["Pos", "EBR"])
    mask_df = mask_df.sort_values(by="Pos")
    mask_df.to_csv("EBR_scores.csv", index=False)


def main():
    read_ambiguous_regions()
    refinded_low_confidence_regions()
    snp_hotspots()
    empiral_base_level_recall()

if __name__ == "__main__":
    main()
