# Biseq mode view in IGV from ONT data

To be launch after stretch_pipeline.sh to have chrX/chrX_methylation.bed file with read informations.


# Known errors

```bash
Sort the file to avoid having position reversed in minus strand with : sort -k1,1 -k4,4 -k2,2n chr1_methylation.bed > chr1_methylation.sort.bed
Error :

###
Chr1 / Chr1_RagTag --> custom name of chr
Custom chr name are not processed correctly for now
###


###
input chrom for all; folder architecture :

../data/chr1/meth_chr1.bed
../data/chr2/meth_chr2.bed
../data/chr3/meth_chr3.bed
../data/chr4/meth_chr4.bed
../data/chr5/meth_chr5.bed

or

../data/chrom_folder/meth_chr1.bed
../data/chrom_folder/meth_chr2.bed
../data/chrom_folder/meth_chr3.bed
../data/chrom_folder/meth_chr4.bed
../data/chrom_folder/meth_chr5.bed
###


###
Polars : 
ImportError : cannot import name 'TypeGuard' 
Solution : pip install typing-extensions --upgrade
###


###
KeyError: '0001e11e-a6bd-4a0b-83b3-989a7bdfe519'
--> Read name not found in fasta file : fastaFromBed didn't work
Solution : Give a genome fasta file with formated lines 

wrong fasta file : 
>Chr1
sequence chr1 all in one line
>Chr2
sequence chr2 all in one line

correct fasta file:
>Chr1
TAAAACCTAAAACCTAAAACCTAAAACCTAAAACCTAAACCCTAAACCCTAAAACCCTAA
ACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAACCTAAACCT 
```
