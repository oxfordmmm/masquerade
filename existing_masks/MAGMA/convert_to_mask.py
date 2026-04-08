import pandas as pd

uvp = pd.read_csv(
    "excluded_loci.txt", sep="\t"
)
print(uvp)
magma = pd.read_csv(
    "UVP_List_of_Excluded_loci.list", header=None, names=["combined"]
)
magma["region"] = magma["combined"].str.split(":", expand=True)[1]
magma["start"] = magma["region"].str.split("-", expand=True)[0].astype(int)
magma["end"] = magma["region"].str.split("-", expand=True)[1].astype(int)
print(magma)

uvp_sites: set[int] = set()
for start, end in zip(uvp["chromStart"], uvp["chromEnd"]):
    # file is 1-indexed inclusive
    # start in inclusive, end is exclusive
    uvp_sites.update(range(start - 1, end))
    
magma_sites: set[int] = set()
for start, end in zip(magma["start"], magma["end"]):
    # file is 1-indexed inclusive
    # start in inclusive, end is exclusive
    magma_sites.update(range(start - 1, end))

with open("uvp.mask", "w", encoding="utf-8") as f:
    for site in sorted(uvp_sites):
        f.write(f"{site}\n")
with open("magma.mask", "w", encoding="utf-8") as f:
    for site in sorted(magma_sites):
        f.write(f"{site}\n")
