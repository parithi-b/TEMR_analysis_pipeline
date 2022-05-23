### vcf files generate using short-read(manta/delly/lumpy) & long-read (pbsv/sniffles/svim) SV callers, and additionally duphold was run to calculate read-depth information.

Tools (version and parameter setting)
#### short-read
<ol>
    <li><a href=https://github.com/lh3/bwa>BWA-MEM</a>(v0.7.17 and default)</li>
    <li><a href=https://github.com/Illumina/manta>Manta</a>(v1.3.2 and High sensitivity calling)</li>
    <li><a href=https://github.com/arq5x/lumpy-sv>LUMPY</a>(v0.2.13 and default)</li>
    <li><a href=https://github.com/dellytools/delly>DELLY</a>(v0.7.8 and default)</li>
    <li><a href=https://github.com/brentp/duphold>Duphold</a>(v0.2.1 and default)</li>
    <li><a href=https://github.com/brentp/mosdepth>mosdepth</a>(v0.3.2 and default)</li>
</ol>

#### long-read
<ol>
    <li><a href=https://github.com/PacificBiosciences/pbh5tools>pbh5tools</a>(v0.8.0 and default)</li>
    <li><a href=https://github.com/philres/ngmlr>NGMLR</a>(v0.2.6 and default)</li>
    <li><a href=https://github.com/fritzsedlazeck/Sniffles>Sniffles</a>(v1.0.7 and default (and -s 5))</li>
    <li><a href=https://github.com/pacificbiosciences/pbsv/>pbsv</a>(v2.2.0 and default)</li>
    <li><a href=https://github.com/eldariont/svim>SVIM</a>(v1.4.0 and default)</li>
    <li><a href=https://github.com/EichlerLab/pav>pav</a>(v1.1.0 and default)</li>
    <li><a href=https://github.com/brentp/duphold>Duphold</a>(v0.2.1 and default)</li>
</ol>
