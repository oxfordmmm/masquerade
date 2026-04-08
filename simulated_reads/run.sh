
# Simulate 150bp reads (though it's a fasta file)
python3 simulate.py ../h37rv.fa --read-length 150 --spacing 10 reads.fasta

# Use blastn
mkdir -p h37rvDB
makeblastdb -in ../h37rv.fa -dbtype nucl -out h37rvDB/db -title "H37Rv DB"

# defaults for blastn task: word size is 11, gapopen is 5, gapextend is 2
# default for megablast task: word size is 28, gapopen is 0, gapextend is 0
# see https://www.ncbi.nlm.nih.gov/books/NBK279684/


blastn -task megablast -query reads.fasta -db h37rvDB/db \
       -dust no -outfmt "6 qseqid qstart qend sseqid sstart send length mismatch pident evalue" \
       -evalue 0.0001 -num_threads 10 \
       -word_size 17 -gapopen 5 -gapextend 2 \
       -max_target_seqs 100000 \
       -out blast.tsv

# Find overlaps

min_length=100
min_identity=90
python3 find_overlaps.py blast.tsv \
    --min-length ${min_length} \
    --min-identity ${min_identity} \
    -m simulated.mask \
    -b filtered_overlaps.tsv