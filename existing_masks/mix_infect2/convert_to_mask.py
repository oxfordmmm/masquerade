import pandas as pd

df = pd.read_csv("MaskedRegions.csv")

mask: set[int] = set()

res_mask: set[int] = set()
pe_ppe_mask: set[int] = set()
repeat_mask: set[int] = set()

for start, stop, cds in zip(df["START"], df["STOP"], df["CDS"]):
    # convert to 0-indexed
    mask.update(range(start - 1, stop))

    if cds == "Drug_resistance":
        res_mask.update(range(start - 1, stop))
    elif cds == "PE/PPE":
        pe_ppe_mask.update(range(start - 1, stop))
    elif cds == "Repeat":
        repeat_mask.update(range(start - 1, stop))

with open("mix_infect2.mask", "w", encoding="utf-8") as f:
    for site in sorted(mask):
        f.write(f"{site}\n")

for file, mask_set in zip(
    ["res.mask", "ppe.mask", "repeat.mask"],
    [res_mask, pe_ppe_mask, repeat_mask],
):
    with open(file, "w", encoding="utf-8") as f:
        for site in sorted(mask_set):
            f.write(f"{site}\n")
