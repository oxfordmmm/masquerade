import pandas as pd

df = pd.read_csv("mhall_mask.bed", header=None, sep="\t", names=["chr", "start", "end"])
df["diff"] = df["end"] - df["start"]
df.to_csv("tmp.tsv", sep="\t", index=False)

sites: set[int] = set()
for start, end in zip(df["start"], df["end"]):
    # start in inclusive, end is exclusive
    sites.update(range(start, end))

sorted_sites = sorted(sites)
with open("mhall.mask", "w", encoding="utf-8") as f:
    for site in sorted_sites:
        f.write(f"{site}\n")
