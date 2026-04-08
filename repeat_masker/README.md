# RepeatMasker

Online tutorials are lacking, but I found that the best use was just using RepeatModeler as follows, following the github repo:
https://github.com/Dfam-consortium/RepeatModeler

You will need to conda install `repeatscout` and `repeatmodeler` to run this.

```bash
mkdir -p h37rv_db
BuildDatabase -name h37rv_db/DB ../h37rv.fa
RepeatModeler -database h37rv_db/DB -threads 40 -LTRStruct -srand 42

RepeatMasker -e ncbi -pa 40 -gff -lib RM_*/consensi.fa.classified -dir MaskerOutput ../h37rv.fa
```
The `output h37rv.fa.out.gff` will have all the masking.

Convert to mask:
```bash
python3 mask_from_gff.py MaskerOutput/h37rv.fa.out.gff --min_score 20
```
Which will produce an overall mask, and one just using the repeat gene families found by RepeatModeler.

---

Repeat masker can also use some database, but it didn't work very well for me.
