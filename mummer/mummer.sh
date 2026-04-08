
nucmer --maxmatch -p mummer_self ../h37rv.fa ../h37rv.fa -t 10

show-coords -rcl -T mummer_self.delta > mummer_self.tsv

exact-tandems ../h37rv.fa 31 > mummer_tandems.txt

python3 mummer_to_mask.py --min_extent 75 --input mummer_self.tsv --tandem_file mummer_tandems.txt mummer.mask
python3 mummer_to_mask.py --min_extent 75 --tandem_file mummer_tandems.txt just_tandems.mask
