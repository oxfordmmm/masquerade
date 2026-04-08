import pandas as pd

df = pd.read_csv(
    "excluding_regions.bed", header=None, sep="\t", names=["chr", "start", "end", "name"]
)

sites: set[int] = set()
for start, end in zip(df["start"], df["end"]):
    # bed file is 0-indexed
    # start in inclusive, end is exclusive
    sites.update(range(start, end))

with open("tbseqpipe.mask", "w", encoding="utf-8") as f:
    for site in sorted(sites):
        f.write(f"{site}\n")
