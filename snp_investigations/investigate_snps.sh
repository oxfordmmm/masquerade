
mhall_mask=/home/ubuntu/studies/masquerade/existing_masks/MHall/mhall.mask
tb_exclude=/home/ubuntu/compare/masks/tb-exclude.txt
dust_mask=/home/ubuntu/studies/masquerade/dust/dust_25.mask
repeat_mask=/home/ubuntu/studies/masquerade/from_gff/repeat_regions.mask
ppe_mask=/home/ubuntu/studies/masquerade/from_gff/ppe_pe.mask
compass_mask=/home/ubuntu/studies/masquerade/existing_masks/Compass/compass.mask
# res_mask=/home/ubuntu/studies/masquerade/existing_masks/mix_infect2/res.mask
blast_overlap_mask=/home/ubuntu/studies/masquerade/blast/compass_settings.mask
mummer=/home/ubuntu/studies/masquerade/mummer/mummer.mask
marin=/home/ubuntu/studies/masquerade/existing_masks/Marin/marin.mask

blast_db=/home/ubuntu/studies/masquerade/blast/h37rvDB/db

# for region in 1637000:1637200
for region in 1637000:1637200 4359164 55552 208320 1779278 1321729 836300:836500 3777542 3501400:3501700 2268720 3281240 4062629
do
    echo "Investigating SNPs in region: $region"
    
    outdir=region_$region/
    mkdir -p $outdir
    python3 investigate_snps.py \
        --range $region \
        --ref_fasta ../h37rv.fa \
        --snps_csv all_fastas/passed_snps.csv \
        --ont_fasta_dir all_fastas/ont \
        --illumina_fasta_dir all_fastas/illumina \
        --output $outdir \
        --blast_db ${blast_db} \
        --masks $mhall_mask $tb_exclude $dust_mask $repeat_mask $ppe_mask $compass_mask $blast_overlap_mask $mummer $marin
done