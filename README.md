This repository includes scripts to integrate short-read (Illumina) and long-read (PacBio) SV calls used in TEMR analysis.

# TEMR analysis pipeline

Transposable element-mediated rearrangements (TEMRs) are a category of structural variants (&ge;50bp) mediated by transposons.

There are three major steps involved in identifying TEMRs from SV calls.

<ol>
  <li>Extract SVs from each caller/algorithm/tool and organize them in tsv format</li>
  <li>Merge SVs from multiple callers</li>
  <ol style="list-style-type: lower-alpha">
    <li>deletions, duplications, and inversions only</li>
    <li>reciprocal overlap using <a href="https://bedtools.readthedocs.io/en/latest/">bedtools</a> (80% <50kbp and 90% for &ge;50kbp)</li>
    <li>parameters used with short-read SV calls:split-reads(SR), paired-read(PR), read-depth(RD)</li>
    <li>parameters used with long-read SV calls:read-support(RS), read-depth</li>    
    <li>read-depth parameters (DHBFC & DHFFC) were calculated using duphold  (check <a href="https://github.com/brentp/duphold">duphold</a> for additional information)</li>
  </ol>
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
<code><ol type="a">
  <li><b>script</b>: step1_temr_sv_vcf_to_tsv.py</li>
  <li><b>input</b>: vcf_file, sample ID , filter/nofilter, short-read/long-read and tool name [manta/delly/lumpy/pbsv/svim/sniffle]</li>
  <li><b>output</b>: tsv files containing SVs calls in the following format</li>
    <ul style="list-style-type: lower-alpha">
      <li>short-read : [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, SR, PR, DHBFC and DHFFC]</li>
      <li>long-read  : [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, RS, DHBFC and DHFFC]</li>
    </ul>
</ol>
###### filtering is applied to remove SVs overlapping simple repeats (50% overlap) and SVs near gaps and centromeres (<500bp)<br>
</code>
  
### STEP 2
<ol>
  <li><b>script</b>: step2_temr_merge_sv_multiple_callers.py</li>
  <li><b>input</b>: </li>
  <li><b>output</b>: </li>
</ol><br>
  
### STEP 3
<ol>
  <li><b>script</b>: step2_temr_merge_sv_multiple_individuals.py</li>
  <li><b>input</b>: </li>
  <li><b>output</b>: </li>
</ol><br>
  
### STEP 4
<ol>
  <li><b>script</b>: step2_temr_identify_transposons_at_breakpoint.py</li>
  <li><b>input</b>: </li>
  <li><b>output</b>: </li>
</ol><br>
