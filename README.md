## Installation

To install geneRFinder, please running the script follow:
``` R
Rscript ./src/config.R
```

## Usage

To run geneRFinder, 
``` R
Rscript ./geneRFinder.R -i [fasta_file_name] -o [output_file_name] -t [thread_number] -s [type_start]
``` 

<br />
[type_start]: <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1 - if start codon is ATG <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2 - if start codon is ATG, GTG and TTG <br />
<br />

For example,

``` R
Rscript ./geneRFinder.R -i ./example/final.contigs.fa -o genes -t 7 -s 1
``` 
