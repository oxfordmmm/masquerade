# Make index
genmap index -F ../h37rv.fa -I index/

# run with different parameters
mkdir -p masks
for e in 2 4; do
    for k in 50 75 100; do
        echo "Processing E: $e, K: $k"
        genmap map -I index -O output_e${e}_k${k} -E $e -K $k --csv
        python3 create_mask.py -i output_e${e}_k${k}.csv -o masks/e${e}_k${k}.txt
    done
done