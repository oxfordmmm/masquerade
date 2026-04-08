# Masquerade - The mask making machine

For our studies evaluating the performance of ONT for sequencing TB we wanted a mask to use when calculating sample relatedness.
Several masks have already been created for slightly different purposes. The main one being the one used by UKHSA in the Compass pipeline (see section "Existing Masks").

However, it's not clear exactly how this mask was created and there are indeed a couple different versions of it.
As such this small project sets out how we created our mask, and uses steps that could be simply applied to other genome references.

**Note**: This mask is only intended for calculating relatedness. It is not used in the assembly stage etc.

The resulting mask is used in the following papers:
- ONT protocol for Mycobacterial MGIT cultures: https://www.biorxiv.org/content/early/2026/02/06/2026.02.04.703726

### Preferred method:
- dust with default parameter (20)
- mummer for tandem repeats
- direct_repeats from ncbi annotation pipeline
- self blast using near compass settings. (megablast, gapopen=5, gapextend=2, min_length=100, min_identity=90)

If using tb specific knowledge:
- marin RLC mask
- repeat_regions from gff file
- list of sites with illumina/ont differences

For any mask it is worth closing up small gaps (<= 100 bp) in the mask.

If being conservative can add:
- RepeatMasker families with score > 20
- genmap using k=50,e=4


### ONT vs Illumina differences

`ont_vs_illumina_snp_density.csv` records all sites with snps between ont and illumina on the same samples (produced as part of the ONT vs Illumina study). This can be used as a way of checking that a generated mask covers those sites of difference.

The best estimate for where Illumina maps incorrectly is from Marin et al (see section "Existing Masks").


## Installation
Everything is installed via the following:

```bash
conda env create -f env.yml
conda activate masquerade
pip install -e .[dev]
```

## Building composite masks
Each method has a directory. May need to run a bash script or python script to make relavant mask before combining.
The table `mask_descriptions.csv` shows how many additional bases each step of filtering adds.

```bash
mkdir -p created_masks

combine_masks dust/dust_20.mask mummer/just_tandems.mask ncbi_annotation_pipeline/ncbi.mask blast/self_blast.mask -o created_masks/auto.mask
# 187130
expand_mask created_masks/auto.mask -e 100 created_masks/auto.mask
# 197876

combine_masks existing_masks/Marin/marin_RLC.mask from_gff/repeat_regions.mask -o created_masks/tb_specific.mask
# 191137

combine_masks created_masks/auto.mask created_masks/tb_specific.mask -o created_masks/masquerade_std_tb.mask
# 260560

expand_mask created_masks/masquerade_std_tb.mask -e 100 created_masks/masquerade_std_tb.mask
# 264525

# manually add remaining problematic sites to produce created_masks/masquerade_std_tb_plus.mask
```

The following script can be used to find the which of the ont_vs_illumina snps are covered by the masks, or how far away they are from the masks. Used to guide mask creation.
```bash
look_up_location --find_distance ont_vs_illumina_snp_density.csv ont_vs_illumina_snp_density_labelled.csv  \
    dust/dust_20.mask mummer/just_tandems.mask ncbi_annotation_pipeline/ncbi.mask blast/self_blast.mask \
    repeat_masker/repeat_families.mask genmap/masks/e4_k50.txt \
    existing_masks/Marin/marin_RLC.mask from_gff/repeat_regions.mask \
    existing_masks/tb-exclude.txt existing_masks/MHall/mhall.mask \
    from_gff/ppe_pe.mask \
    created_masks/masquerade_std_tb.mask
```

## Note on file formats

Bed files are 0-indexed with start being inclusive, and end being exclusive. (https://en.wikipedia.org/wiki/BED_(file_format))

GFF3 are 1-indexed with both start and stop being inclusive

Generated mask files are just lists of sites and are 0-indexed

---

## Existing masks
This is an incomplete summary of the TB masks I've found.

### TB-Exclude FN5 (0-indexed)
This comes from FN4 which is supposedly based on compass, but unknown origins aobut how it was created.
It is 0-indexed.

### Compass Mask (1-indexed into 0-indexed)
Tim Peto provided a file which is meant to be the compass mask, and presumably comes from [compass](https://github.com/oxfordmmm/CompassCompact/tree/master), but not sure of the link. This looks 1-indexed.

The mask has separately marked the rrl and rrs sections.

I've converted this to `existing_masks/compass_mask.txt` which is 0-indexed.

I've been unable to run compass, but have tried running the blastn command using the provided settings but get a much smaller mask than that provided by Tim.

### Michael Hall (Bed file into 0-indexed list)
In this [paper](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(22)00301-9/fulltext) he looks at clustering with ONT vs Illumina.

He uses the mask `existing_masks/mhall_mask.bed` from [here](https://github.com/mbhall88/head_to_head_pipeline/blob/master/analysis/baseline_variants/resources/compass-mask.bed)

He says that he got the mask from a [Tim walker paper](https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(14)70027-X/fulltext) but I can't find it in here.


### Comparing

The mask from Tim Peto is strictly larger than the Michael Hall mask (by 7293), mainly it seems due to the rrl (3138) and rrs (1537) sites.

The tb-exclude is missing 49567 sites which are present in the Peto mask, but has 5372 additional sites

Combining all together covers 337636 sites.


### Marin - based on emprical accuracy of short reads
[paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8963317)

Is only 177kbp long and is built based on where Illumina has poor mapping.

Files have been copied from supplementary data from the paper.

Marin has several systems for masking:
- EBR (emprial base-level recall) is proportion of isolates making a confident+correct variant call at site.
- pileup mappability (using genmap at k=50, e=4) for 50-mers which matches (allowing for 4 errors)
- PLC is putative low confidence. PPE/PE genes, mobile genetic elements and 69 additional genes with homology elsewhere. (10% of genome)

They produce a set of refined low confidence regions (RLC) - S13. Made from:
- false snp hot spots (S9/21)
- EBR < 0.9 (S18)
- regions with ambigously defined by long read sequencing (S5)



### Modlin - blind spots in illumina sequence 
[paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8190613/#s5)
Idea is that they found low coverage parts of genome which vary by illumina platform used.
Could be used as a kind of test set for illumina side of things.

Requires the `S2-pooled-blindspots.csv` from the supplement.


### MixInfect2
This [repo](https://github.com/bensobkowiak/MixInfect2?tab=readme-ov-file) is for detecting mixed strain infections.
It contains a mask file which I think it is coming directly from the gff3 file for h37rv, which labels genes, PPE and repeat regions.
It does helpfully label the genes/reason for masking a region.

By comparison to gff3 I think the start and stop are 1-indexed inclusive.

Length of mask is 362535, which breaks down to:
- Res genes: 54021
- PPE/PE: 286701
- Repeat regions: 23862

Notes:
- They have included ispE is the PE family (I think a mistake)
- wag22/esxS/lipY is actually part of PE family

### TBSeqPipe
From this [repo](https://github.com/KevinLYW366/TBSeqPipe/tree/master/).
Doesn't seem to have a detailed paper.
The bed file labels the genes it is using, which mostly maps to the gff3 file.
But names used are the gene ids rather than the more descriptive ones (e.g. `Rv0755c`).

Mask is very large at 482234.

Does mask mobile genetic elements which might be smart.


### MAGMA
[MAGMA repo](https://github.com/TORCH-Consortium/MAGMA) uses mask from previous pipeline [UVP](https://github.com/CPTR-ReSeqTB/UVP).

MAGMA uses [UVP_List_of_Excluded_loci.list](https://github.com/TORCH-Consortium/MAGMA/blob/master/resources/regions/UVP_List_of_Excluded_loci.list.)
I think the original file from UVP comes from [here](https://github.com/CPTR-ReSeqTB/UVP/blob/master/uvp/data/excluded_loci.txt).

The UVP file labels the regions, but MAGMA doesn't.
MAGMA also seems to have 4 extra regions:
- 850341-850526
- 1999141-1999356
- 3378328-3378414
- 3380682-3380992
Which don't cleanly match the gff3 file.

Magma mask is 462123bp.

The indexing matches the gff3 file. So 1-indexed inclusive regions.

### MTBseq
[repo](https://github.com/ngs-fzb/MTBseq_source).

Not clear how they are masking.
They do have a [file](https://github.com/ngs-fzb/MTBseq_source/blob/master/var/cat/MTB_Gene_Categories.txt) which labels genes as being essential, nonessential or repetitive. 
