This repository includes scripts to integrate short-read (Illumina) and long-read (PacBio) SV calls used in TEMR analysis.

# TEMR analysis pipeline

Transposable element-mediated rearrangements (TEMRs) are a category of structural variants (&ge;50bp) mediated by transposons.

There are three major steps involved in identifying TEMRs from SV calls.
<ol>
  <li>Extract SVs from each caller/algorithm/tool and organize them in tsv format</li>
  <li>Merge calls from multiple callers</li>
    <dl>SV types currently considered: deletions, duplications, and inversion</dl>
    <dl>SV type future consideration: insertion</dl>
  <li>Merge calls from multiple individuals</li>
  <li>Identify SVs mediated by TEs</li>
</ol>
