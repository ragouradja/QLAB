# Biseq mode view in IGV from ONT data

To be launch after stretch_pipeline.sh to have chrX/chrX_methylation.bed file with read informations.


## Slow and Fast mode
### Slow mode
Mode by default. Will read all 5 chromosomes files in parallel (if `--all` is used) without loading the data into memory. We save RAM memory (using around 0%) but cost more time.

### Fast mode
Need to be activated with `--fast`. Wille read all 5 chromosomes files in parallel (if `--all` is used) and will load the data into memory. We save time but cost RAM memory (at least 3x the total size of `chrX/chrX_methylation.bed`)

#### Memory usage in fast mode
<p align="center">
<img src="https://user-images.githubusercontent.com/71189947/190190992-d51145f1-59d9-4d94-91b9-e92708fe2214.png" width="500" height="500"/>
</p>

### Time comparaison
<p align="center">
<img src="https://user-images.githubusercontent.com/71189947/190188112-dae6ee8a-211f-4118-a64f-2f35eb4d692e.png" width="500" height="500"/>
</p>

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
