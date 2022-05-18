This repository includes scripts to integrate short-read (Illumina) and long-read (PacBio) SV calls used in TEMR analysis.

# TEMR analysis pipeline

Transposable element-mediated rearrangements (TEMRs) are a category of structural variants (&ge;50bp) mediated by transposons.

There are four major steps involved in identifying TEMRs from SV calls.

<ol>
  <li>Extract SVs from each caller/algorithm/tool and organize them in tsv format</li>
   <ol style="list-style-type: lower-alpha">
    <li>deletions, duplications, and inversions only</li>
    <li>additional information for short-read SV calls:  paired-read(PR), split-reads(SR), read-depth(RD)</li>
    <li>additional information for long-read SV calls: read-support(RS), read-depth</li>    
    <li>read-depth parameters (DHBFC & DHFFC) were calculated using duphold  (check <a href="https://github.com/brentp/duphold">duphold</a> for additional information)</li>
  </ol>
  <li>Merge SVs from multiple callers</li>
  <ol style="list-style-type: lower-alpha">
    <li>reciprocal overlap using <a href="https://bedtools.readthedocs.io/en/latest/">bedtools</a> (80% <50kbp and 90% for &ge;50kbp)</li>   
  <li>Merge SVs from multiple individuals</li>
  <ol style="list-style-type: lower-alpha">
    <li>reciprocal overlap using <a href="https://bedtools.readthedocs.io/en/latest/">bedtools</a> (80% <50kbp and 90% for &ge;50kbp)</li>  </ol>    
  <li>Identify SVs mediated by TEs</li>
  <ol style="list-style-type: lower-alpha">
    <li>RepeatMasker from UCSC table browser</li>
  </ol>
</ol>

###### please use the full path of the files and have all the required scripts & vcfs in a single folder for ease of use. Python and bedtools were run using conda environment. 

### STEP 1
<ol type="a">
  <li><b>script</b>: step1_temr_sv_vcf_to_tsv.py</li>
  <li><b>input</b>: vcf_file, sample ID , filter/nofilter, short-read/long-read and tool name [manta/delly/lumpy/pbsv/svim/sniffle]</li>
  <li><b>output</b>: tsv files containing SVs calls in the following format</li>
    <ul style="list-style-type: lower-alpha">
      <li>short-read : [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, PR, SR, DHBFC and DHFFC]</li>
      <li>long-read  : [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, RS, DHBFC and DHFFC]</li>
    </ul>
</ol>

###### filtering is applied to remove SVs overlapping simple repeats (50% overlap) and SVs near gaps and centromeres (<500bp).

```
Example 
input:
python step1_temr_sv_vcf_to_tsv.py vcf_files/HG00733_manta_duphold.vcf HG00733 filter short-read manta

output: 
filename: vcf_files/NA19240_manta_duphold_sv_filtered.tsv
(sample SV from the output file)
chr1	9226573	9228034	DEL	NA19240	manta	36	22	0.40625	0.382353
```
  
### STEP 2
<ol>
  <li><b>script</b>: step2_temr_merge_sv_multiple_callers.py</li>
  <li><b>input</b>: three tsv files, corresponding caller names(in order on the files mentioned), parameters (PR,SR,RD or RS,RD) </li>
  <li><b>output</b>: </li> tsv file containing merged SV calls
      <ul style="list-style-type: lower-alpha">
      <li>short-read : [CHR, POS, END, SVTYPE, SAMPLEID, CALLERS]</li>
    </ul>
</ol><br>
  
```
Example 
input:
python step2_temr_merge_sv_multiple_callers.py vcf_files/NA19240_manta_duphold.vcf
vcf_files/NA19240_delly_duphold.vcf vcf_files/NA19240_lumpy_duphold.vcf NA19240 10 5 True

If an SV fails to reach the required support it is filtered out before merging
10 -> minimum number of paired-read support needed
5 -> minimum number of split-read support needed
True -> default RD value from duphold

few example combination of parameters for filtering
[10 5 False] --> only PR and SR filter
[0 0 True] --> only RD filter
[0 0 False] --> default (no filter)
output: 
filename: vcf_files/Ensemble/NA19240/NA19240_10_5_RD_merged_sorted_under50k.bed
(sample SV from the output file)
chr1	9226573	9228034	DEL	NA19240	manta;delly;lumpy
```

###### We merge SVs that passed the filtering step using bedtools. The default value is set to 80% RO (used in this study), but can be changed if needed. Additionally, SVs from multiple callers are merged using a Rank based method [1. manta, 2.delly, 3. lumpy  -> was used in this study]. The SV from a highest ranked caller is retained while the rest are excluded. For example, if an SV is called by manta, delly and lumpy, during the merging process manta call is retained while delly and lumpy calls are removed. [More information can be found in the manuscript]

### STEP 3
<ol>
  <li><b>script</b>: .py</li>
  <li><b>input</b>: </li>
  <li><b>output</b>: </li>
</ol><br>
  
### STEP 4
<ol>
  <li><b>script</b>: .py</li>
  <li><b>input</b>: </li>
  <li><b>output</b>: </li>
</ol><br>
