This repository includes scripts to integrate short-read (Illumina) and long-read (PacBio) SV calls used in TEMR analysis.

# TEMR analysis pipeline

Transposable element-mediated rearrangements (TEMRs) are a category of structural variants (&ge;50bp) mediated by transposons.

There are three major steps involved in identifying TEMRs from SV calls.

<ol>
  <li>Extract SVs from each caller/algorithm/tool and organize them in tsv format</li>
  <li>Merge SVs from multiple callers</li>
  <ol style="list-style-type: lower-alpha">
    <li>deletions, duplications, and inversions only</li>
    <li>reciprocal overlap (80% <50kbp and 90% for &ge;50kbp)</li>
    <li>parameters used with short-read SV calls:split-reads(SR), paired-read(PR), read-depth(RD)</li>
    <li>parameters used with long-read SV calls:read-support(RS), read-depth</li>    
    <li>read-depth are calcualted using duphold DHBFC & DHFFC (check https://github.com/brentp/duphold)</li>

  </ol>
  <li>Merge SVs from multiple individuals</li>
  <ol style="list-style-type: lower-alpha">
    <li>reciprocal overlap (80% <50kbp and 90% for &ge;50kbp)</li>
  <li>Identify SVs mediated by TEs</li>
</ol>

###### please use the full path of the files and please all the required scripts & vcfs in a single folder for ease of use.
  
####STEP 1
script : step1_temr_vcf_to_be.py
input : vcf_file, sample ID , short-read/long-read and tool name [manta/delly/lumpy/pbsv/svim/sniffle]
ooutput : tsv files containing SVs calls in the following format 
  short-read: [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, SR, PR, DHBFC and DHFFC]
  long-read: [CHR, POS, END, SVTYPE, SAMPLEID, CALLER, RS, DHBFC and DHFFC]
