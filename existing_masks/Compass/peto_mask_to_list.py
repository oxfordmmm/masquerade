import pandas as pd

df = pd.read_csv("compass_mask_from_peto.csv")

print(df["gene"].unique())

df = df[df["gene"].isin(["mask", "rrs", "rrl"])]
print(df)

# convert to 0-indexed
df["pos"] = df["pos"].astype(int) - 1

df = df[["pos"]]
df = df.sort_values(by="pos")
df.to_csv("compass.mask", index=False, header=False)