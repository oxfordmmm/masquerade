

mkdir -p overlaps
for l in 75 100 125; do
    for i in 80 85 90 95; do
        echo "Running find_overlaps.py with min-length $l and min-identity $i"
        python3 find_overlaps.py blast.tsv \
            --min-length $l \
            --min-identity $i \
            -m overlaps/${l}_${i}.mask \
            -b overlaps/${l}_${i}.tsv
    done
done