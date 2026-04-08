mkdir -p h37rvDB
makeblastdb -in ../h37rv.fa -dbtype nucl -out h37rvDB/db -title "H37Rv DB"

# defaults for blastn task: word size is 11, gapopen is 5, gapextend is 2
# default for megablast task: word size is 28, gapopen is 0, gapextend is 0
# see https://www.ncbi.nlm.nih.gov/books/NBK279684/

# compass compact use megablast with gapopen 5, gapextend 2, word size 17

# reducing word size increases sensitivity but also runtime, so using 15


# Uses the settings from compass compact
blastn -task megablast -query ../h37rv.fa -db h37rvDB/db \
       -dust no -outfmt "6 qseqid qstart qend sseqid sstart send length mismatch pident evalue" \
       -evalue 0.0001 -num_threads 20 \
       -word_size 15 -gapopen 5 -gapextend 2 \
       -max_target_seqs 10000000 \
       -out self_blast.tsv


# Compass mask used length 75 and identity 90
# Have increased min length to 100
python3 blast_to_mask.py self_blast.tsv self_blast.mask --min_length 100 --min_identity 90



# WARNING: This does not seem to reproduce the original compass mask!