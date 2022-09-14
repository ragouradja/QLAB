# Biseq mode view in IGV from ONT data

Creating BAM file with Bisulfite treatment simulated from ONT BED file.

To be launch after stretch_pipeline.sh to have chrX/chrX_methylation.bed file with read informations.

## Usage
```bash
$ python ont_to_bam.py --help

usage: ont_to_bam.py [-h] --file chr_file [--genome genome] --basename_output BASENAME_OUTPUT [--output OUTPUT] [--fast] [--all]

Generate bam file from ONT reads data

optional arguments:
  -h, --help            show this help message and exit
  --file chr_file, -f chr_file
                        Path to file with reads informations for one chromosome. (default: None)
  --genome genome, -g genome
                        Path to the genome fasta file (Default : Col-CEN) (default: /mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta)
  --basename_output BASENAME_OUTPUT, -base BASENAME_OUTPUT
                        Output file without extension (default: None)
  --output OUTPUT, -o OUTPUT
                        Path to output directory. (default: .)
  --fast                Loading data into memory to work faster (Relatively high memory usage, need at least 2.5*FILE_SIZE Go of RAM). (default: False)
  --all                 Compute all chromosomes in parallel. (default: False)
```

By default, the script will run in slow mode for one chromosome given in `--file` and will use Col-CEN genome (you can change it with `--genome`).


```bash
$ python ont_to_bam.py --file 
```


## Slow and Fast mode
### Slow mode
Mode by default. Will read all 5 chromosomes files in parallel (if `--all` is used) without loading the data into memory. We save RAM memory (using around 0%) but cost more time.

### Fast mode
Need to be activated with `--fast`. Wille read all 5 chromosomes files in parallel (if `--all` is used) and will load the data into memory. We save time but cost RAM memory (around 2.5x the total size of BED files)

#### Memory usage in fast mode
<p align="center">
<img src="https://user-images.githubusercontent.com/71189947/190192041-2315610f-ffa7-4cf9-b6e3-6e2b5b1d3623.png" width="700" height="500"/>
</p>

### Time comparaison
<p align="center">
<img src="https://user-images.githubusercontent.com/71189947/190188112-dae6ee8a-211f-4118-a64f-2f35eb4d692e.png" width="700" height="500"/>
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

ModuleNotFoundError: No module named 'polars'
Solution :  pip install polars

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
