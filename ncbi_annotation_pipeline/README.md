# NCBI annotation pipeline

This follows steps from https://github.com/ncbi/pgap/wiki/Quick-Start.

```bash
wget -O pgap.py https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py
chmod +x pgap.py
./pgap.py --update

./pgap.py -r -o results -g ../h37rv.fa -s 'Mycobacterium tuberculosis'
```

Then extract the `direct_repeat` sections using
```bash
python3 mask_from_gff.py results/annot.gff ncbi.mask
```
