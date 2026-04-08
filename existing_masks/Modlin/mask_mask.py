import pandas as pd

df = pd.read_csv("S2-pooled-blindspots.csv")

df = df[["blindspot"]]
# make 0-indexed
df["blindspot"] = df["blindspot"] - 1
df.to_csv("blindspots.mask", index=False, header=False)
