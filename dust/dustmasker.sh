
# level 20 is the default for dustmasker

level=$1
if [ -z "$level" ]; then
  echo "Usage: $0 <dust level>"
  exit 1
fi

dustmasker -in ../h37rv.fa -infmt fasta -parse_seqids \
  -outfmt interval -out dust_${level}.tsv -level ${level}

python3 dust_to_mask.py dust_${level}.tsv dust_${level}.mask

rm dust_${level}.tsv

# window masker seems to mask to much??
# windowmasker -in h37rv.fa -infmt fasta -parse_seqids \
#   -mk_counts -out dust_counts.tsv
# windowmasker -in h37rv.fa -infmt fasta -parse_seqids \
#   -ustat dust_counts.tsv -dust T -outfmt interval -out dust_window.tsv
