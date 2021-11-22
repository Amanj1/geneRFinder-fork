## Installation

To install geneRFinder, please running the script follow:
``` R
Rscript ./src/config.R
```

## Usage

To run geneRFinder, 
``` R
Rscript ./geneRFinder.R -i [fasta_file_name] -o [output_file_name] -t [thread_number] -s [start_type] -n [intergenic]
``` 

<br />
[fasta_file_name]: input file name

<br />
[output_file_name]: output file name

<br />
[thread_number]: number of thread

<br />
[start_type]: <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1 - if start codon is ATG <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2 - if start codon is ATG, GTG and TTG <br />

<br />
[intergenic]: <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1 - output without intergenic sequences <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2 - output with intergenic sequences <br />
<br />

For example,

``` R
Rscript ./geneRFinder.R -i ./example/final.contigs.fa -o output -t 7 -s 1 -n 1
``` 

Please, download the src/model.RData file separately, it is a large file.


## Usage

In the output there are for each sequence:
```
>[random_id], [sequence_length],[contig_id] <br />
[sequence]
```

Example:
```
>34, len=63, k121_33 flag=1 multi=2.0000 len=501
ATGATAAAAGCGCGCGTCAGGTACGGCTCGTCGCCGCCGGCAATGCCTATGCGGTCACGCTAA
```
