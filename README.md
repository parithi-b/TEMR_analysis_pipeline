This repository includes scripts to integrate short-read (Illumina) and long-read (PacBio) SV calls used in TEMR analysis.

# TEMR analysis pipeline

Transposable element-mediated rearrangements (TEMRs) are a category of structural variants (&ge;50bp) mediated by transposons.

There are five major steps involved in identifying TEMRs from SV calls.

<ol>
  <li>Extract SVs data from each caller/algorithm/tool and organize them in tsv format</li>
   <ol style="list-style-type: lower-alpha">
    <li>deletions, duplications, and inversions only</li>
    <li>additional information for short-read SV calls:  paired-read(PR), split-reads(SR), read-depth(RD)</li>
    <li>additional information for long-read SV calls: read-support(RS), read-depth</li>    
    <li>read-depth parameters (DHBFC & DHFFC) were calculated using duphold  (check <a href="https://github.com/brentp/duphold">duphold</a> for additional information)</li>
  </ol>
  
  <li>Filter and Merge SVs from multiple callers </li>
  <ol style="list-style-type: lower-alpha">
    <li>reciprocal overlap (80%) using <a href="https://bedtools.readthedocs.io/en/latest/">bedtools</a> </li>   
  </ol> 
  
  <li>Merge SVs from multiple individuals</li>
  <ol style="list-style-type: lower-alpha">
    <li>reciprocal overlap (80%) using <a href="https://bedtools.readthedocs.io/en/latest/">bedtools</a> </li>  </ol>    
    
  <li>Merge SVs from long-read ensemble and short-read ensemble</li>
  <ol style="list-style-type: lower-alpha">
    <li>reciprocal overlap (80%) using <a href="https://bedtools.readthedocs.io/en/latest/">bedtools</a> </li>  </ol> 
    
  <li>Identify SVs mediated by TEs</li>
  <ol style="list-style-type: lower-alpha">
    <li>RepeatMasker from UCSC table browser</li>
  </ol>
</ol>

###### please use the full path of the files and have all the required scripts & vcfs in a single folder for ease of use. Python and bedtools were run using conda environment. 

---

### STEP 1: Extract SV data 

###### In this step we filter SVs overlapping simple repeats (50% overlap) and SVs near gaps and centromeres (<500bp).

<ol type="a">
  <li><b>script</b>: step1_temr_sv_vcf_to_tsv.py</li>
  <li><b>input</b>: vcf_file, sample ID , filter/nofilter, short-read/long-read and tool name [manta/delly/lumpy/pbsv/svim/sniffle]</li>
  <li><b>output</b>: tsv files containing SVs calls in the following format</li>
    <ul style="list-style-type: lower-alpha">
      <li>short-read : [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, PR, SR, DHBFC and DHFFC]</li>
      <li>long-read  : [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, RS, DHBFC and DHFFC]</li>
    </ul>
</ol>

```
Example 
input:
python scripts/step1_temr_sv_vcf_to_tsv.py vcf_files/short-read/NA19240_manta_duphold.vcf NA19240 filter short-read manta

output: 
filename: vcf_files/NA19240_manta_duphold_sv_filtered.tsv
(sample SV from the output file)
chr10	26710086	26713224	DEL	NA19240	manta	127	55	0	0

Same SV identified by other callers
<short-read pipeline>
chr10	26710086	26713225	DEL	NA19240	delly	126	10	0	0
chr10	26710131	26713193	DEL	NA19240	lumpy	125	0	0	0
<long-read pipeline>
chr10	26710086	26713224	DEL	NA19240	pbsv	19	0	0
chr10	26710086	26713224	DEL	NA19240	svim	24	0	0
chr10	26710161	26713284	DEL	NA19240	sniffles	29	0	0
```
  
### STEP 2: Filter and Merge SVs from multiple callers

###### In this step we merge SVs that passed the filters, using bedtools. The default value is set to 80% RO (used in this study), but can be changed if needed. 

###### Additionally, SVs from multiple callers (at least 2 caller support) were merged using a Rank based method [1. manta, 2.delly, 3. lumpy  -> was used in this study]. The SV call from a highest ranked caller is retained while the rest are excluded. For example, if an SV is called by manta, delly and lumpy, during the merging process manta call is retained while delly and lumpy calls are removed. [More information can be found in the manuscript]

<ol>
  <li><b>script</b>: step2_temr_filter_merge_sv_multiple_callers.py</li>
  <li><b>input</b>: three tsv files, corresponding caller names(in order on the files mentioned), parameters (PR,SR,RD or RS,RD), and output folder</li>
  <li><b>output</b>: tsv file containing merged SV calls </li>
      <ul style="list-style-type: lower-alpha">
      <li>[CHR, POS, END, SVTYPE, SAMPLEID, CALLERS]</li>
    </ul>
</ol><br>
  
```
Example 
input:
python step2_temr_merge_sv_multiple_callers.py vcf_files/NA19240_manta_duphold.vcf vcf_files/NA19240_delly_duphold.vcf vcf_files/NA19240_lumpy_duphold.vcf NA19240 10 5 True vcf_files/short_read/Ensemble/

If an SV fails to reach the required support it is filtered out before merging
10 -> minimum number of paired-read support needed
5 -> minimum number of split-read support needed
True -> default RD value from duphold

few example combination of parameters for filtering
[10 5 False] --> only PR and SR filter
[0 0 True] --> only RD filter
[0 0 False] --> default (no filter)
output: 
filename: vcf_files/short_read/Ensemble/NA19240/NA19240_10_5_RD_merged_sorted_under50k.bed
(sample SV from the output file)
chr10	26710086	26713224	DEL	NA19240	manta;delly;lumpy

Same SV in <long-read pipeline>
chr10	26710086	26713224	DEL	NA19240	pbsv;sniffles;svim
```

### STEP 3: Merge SVs from multiple individuals
###### In this step we merge SV calls from 3 individuals
<ol>
  <li><b>script</b>: step3_temr_merge_individuals.py</li>
  <li><b>input</b>: merged SV callset (multicaller merge) from the three individuals, and output folder </li>
  <li><b>output</b>: tsv file containing merged SV calls </li>
      <ul style="list-style-type: lower-alpha">
      <li>[CHR, POS, END, SVTYPE, INDIVIDUALS]</li>
</ol><br>

```
Example
input: 
python step3_temr_merge_individuals.py vcf_files/short_read/Ensemble/HG00514/HG00514_10_5_RD_merged_sorted_under50k.bed vcf_files/short_read/Ensemble/HG00733/HG00733_10_5_RD_merged_sorted_under50k.bed vcf_files/short_read/Ensemble/NA19240/NA19240_10_5_RD_merged_sorted_under50k.bed
vcf_files/short_read/

output:
filename: vcf_files/short_read/All_samples_merged.tsv
(sample SV from the output file)
chr10	26710086	26713224	DEL	HG00514:HG00514;HG00733;NA19240
HG00514: --> lead sample containing this SV
HG00514;HG00733;NA19240 --> all samples contianing this SV

Same SV in <long-read pipeline>
chr10	26710086	26713224	DEL	HG00514:HG00514;HG00733;NA19240
```
  
### STEP 4: Merge SVs from multiple technology 
###### In this step we merge SV calls from both short-read ensemble pipeline and long-read ensemble pipeline. Ensemble pipeline --> multiple caler and multiple individual merge 

<ol>
  <li><b>script</b>: step3_temr_merge_technology.py</li>
 <li><b>input</b>: Ensemble SV callset from short-read pipleline and long-read pipeline, and output folder </li>
  <li><b>output</b>: tsv file containing merged SV calls </li>
      <ul style="list-style-type: lower-alpha">
      <li>[CHR, POS, END, SVTYPE, INDIVIDUALS, SHARED]</li>
</ol><br>

```
Example
input: python step3_temr_merge_technology.py vcf_files/short_read/All_samples_merged.tsv vcf_files/long_read/All_samples_merged.tsv vcf_files
output:
filename: vcf_files/All_samples_short_long_merged.tsv
chr10	26710086	26713224	DEL	HG00514:HG00514;HG00733;NA19240 shared
shared -> identified by both long-read and short-read callers
```

### STEP 5: Identify TEMRs
###### In this step we identify SVs with breakpoints present within two distinct TEs. There is an option to assign a window size while searching for TE, if breakpoint accuracy was a concern.

<ol>
  <li><b>script</b>: step5_temr_identify_te_at_junction.py</li>
  <li><b>input</b>: Final merged file, repeatMasker file</li>
  <li><b>output</b>: tsv file containing TEs present at the breakpoint junction of SV calls </li>
      <ul style="list-style-type: lower-alpha">
      <li>[CHR, POS, END, SVTYPE, INDIVIDUALS, SHARED, 3_TE, 5_TE, TEMR_FLAG]</li>
</ol><br>

```
Example
input:
  python vcf_files/All_samples_short_long_merged.tsv UCSC_track/hg38_repeatMasker.tsv
output:
  filename: vcf_files All_samples_short_long_merged_TEMR.stv
chr10	26710086	26713224	DEL	HG00514:HG00514;HG00733;NA19240	shared	chr10;26709870;26710166;AluSx3;SINE;+;Alu	chr10;26713011;26713314;AluSc8;SINE;+;Alu	TEMR_SAME;Alu
```
