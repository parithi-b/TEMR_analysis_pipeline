## UCSC tracks used in this pipeline
###### please unzip all the files after downloading
<pre><code>
merge all TE files into a single file using the following query
cat *repeatMasker* | sort -k1,1 -k2,2n > hg38_repeatMasker_TEs.tsv
</code></pre>
