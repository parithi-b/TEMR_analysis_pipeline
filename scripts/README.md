## Scripts for running TEMR pipeline

After activating the conda environment the following queries will produce the desired results. (These queires can always be copy to a single bash script and run from there)

---

## Queries 

please run them from the main folder (or adjust the file paths accordingly)

### vcf to tsv
#### short-read 

```
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/NA19240_manta_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv NA19240 manta 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/NA19240_delly_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv NA19240 delly 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/NA19240_lumpy_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv NA19240 lumpy 50 50000

python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/HG00514_manta_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00514 manta 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/HG00514_delly_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00514 delly 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/HG00514_lumpy_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00514 lumpy 50 50000

python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/HG00733_manta_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00733 manta 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/HG00733_delly_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00733 delly 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/HG00733_lumpy_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00733 lumpy 50 50000
```

#### long-read

```
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/NA19240_pbsv_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv NA19240 pbsv 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/NA19240_sniffles_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv NA19240 sniffles 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/NA19240_svim_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv NA19240 svim 50 50000

python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/HG00514_pbsv_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00514 pbsv 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/HG00514_sniffles_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00514 sniffles 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/HG00514_svim_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00514 svim 50 50000

python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/HG00733_pbsv_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00733 pbsv 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/HG00733_sniffles_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00733 sniffles 50 50000
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/long-read/HG00733_svim_duphold.vcf UCSC_tracks/hg38_gaps_centromeres.tsv UCSC_tracks/hg38_simpleRepeats_collapsed.tsv HG00733 svim 50 50000
```

---

### merge multiple caller per individual

#### short-read 

```
python scripts/step2_temr_filter_merge_sv_multiple_callers.py vcf_files/short-read/NA19240_manta_duphold_sv_filtered.tsv vcf_files/short-read/NA19240_delly_duphold_sv_filtered.tsv vcf_files/short-read/NA19240_lumpy_duphold_sv_filtered.tsv manta delly lumpy NA19240 short-read 10 5 True vcf_files/short-read
python scripts/step2_temr_filter_merge_sv_multiple_callers.py vcf_files/short-read/HG00733_manta_duphold_sv_filtered.tsv vcf_files/short-read/HG00733_delly_duphold_sv_filtered.tsv vcf_files/short-read/HG00733_lumpy_duphold_sv_filtered.tsv manta delly lumpy HG00733 short-read 10 5 True vcf_files/short-read
python scripts/step2_temr_filter_merge_sv_multiple_callers.py vcf_files/short-read/HG00514_manta_duphold_sv_filtered.tsv vcf_files/short-read/HG00514_delly_duphold_sv_filtered.tsv vcf_files/short-read/HG00514_lumpy_duphold_sv_filtered.tsv manta delly lumpy HG00514 short-read 10 5 True vcf_files/short-read
```

#### long-read 

```
python scripts/step2_temr_filter_merge_sv_multiple_callers.py vcf_files/long-read/NA19240_pbsv_duphold_sv_filtered.tsv vcf_files/long-read/NA19240_sniffles_duphold_sv_filtered.tsv vcf_files/long-read/NA19240_svim_duphold_sv_filtered.tsv pbsv sniffles svim NA19240 long-read 5 True vcf_files/long-read
python scripts/step2_temr_filter_merge_sv_multiple_callers.py vcf_files/long-read/HG00733_pbsv_duphold_sv_filtered.tsv vcf_files/long-read/HG00733_sniffles_duphold_sv_filtered.tsv vcf_files/long-read/HG00733_svim_duphold_sv_filtered.tsv pbsv sniffles svim HG00733 long-read 5 True vcf_files/long-read
python scripts/step2_temr_filter_merge_sv_multiple_callers.py vcf_files/long-read/HG00514_pbsv_duphold_sv_filtered.tsv vcf_files/long-read/HG00514_sniffles_duphold_sv_filtered.tsv vcf_files/long-read/HG00514_svim_duphold_sv_filtered.tsv pbsv sniffles svim HG00514 long-read 5 True vcf_files/long-read
```
---
### merge multiple individuals per technology 

#### short-read
```
python scripts/step3_temr_merge_individuals.py vcf_files/short-read/Ensemble/NA19240/NA19240_10_5_RD_merged_sorted.tsv vcf_files/short-read/Ensemble/HG00733/HG00733_10_5_RD_merged_sorted.tsv vcf_files/short-read/Ensemble/HG00514/HG00514_10_5_RD_merged_sorted.tsv vcf_files/short-read/
```
#### long-read
```
python scripts/step3_temr_merge_individuals.py vcf_files/long-read/Ensemble/NA19240/NA19240_default_RD_merged_sorted.tsv vcf_files/long-read/Ensemble/HG00733/HG00733_default_RD_merged_sorted.tsv vcf_files/long-read/Ensemble/HG00514/HG00514_default_RD_merged_sorted.tsv vcf_files/long-read/
```
---
### merge two techonologies 
```

```
---
### identify TEMRs
```
```



