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
    <li>parameters used with short-read SV calls:split-reads, paired-read, read-depth</li>
    <li>parameters used with long-read SV calls:read-support, read-depth</li>
  </ol>
  <li>Merge SVs from multiple individuals</li>
  <ol style="list-style-type: lower-alpha">
    <li>reciprocal overlap (80% <50kbp and 90% for &ge;50kbp)</li>
  <li>Identify SVs mediated by TEs</li>
</ol>
