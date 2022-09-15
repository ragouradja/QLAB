# Creating BAM file from ONT BED file to visualize methylation state in IGV.


**Run after stretch_pipeline.sh to get the file chrX/chrX_methylation.bed with the reading information.**

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

Generate BAM file for chromosome 1 in slow mode :
```bash
$ python ont_to_bam.py  --file chr1/chr1_methylation.bed -base col0_test
```

Generate BAM file for all chromosomes if by changing 'chr1' by 'chrX' (X being 2,3,4,5) we can reach the BED file of all 5 chromosomes :
```bash
$ python ont_to_bam.py  --file chr1/chr1_methylation.bed -base col0_test --all
```

Generate BAM file for all chromosomes in Fast mode :

```bash
$ python ont_to_bam.py  --file chr1/chr1_methylation.bed -base col0_test --all --fast
```

Generate BAM file for all chromosomes in Fast mode :
```bash
$ python ont_to_bam.py  --file chr1/chr1_methylation.bed -base col0_test --all --fast
```

If you have all chromosomes in one file (with small data), use `--all_in_one` instead of `--all` :
```bash
$ python ont_to_bam.py  --file col_test.bed -base all_one --all_in_one 
```

## Slow and Fast mode
### Slow mode
Mode by default. Will read all 5 chromosomes files in parallel (if `--all` is used) without loading the data into memory. We save RAM memory (using around 0%) but cost more time.

### Fast mode
Need to be activated with `--fast`. Will read all 5 chromosomes files in parallel (if `--all` is used) and will load the data into memory. We save time but cost RAM memory (around 2.5x the total size of BED files)

#### Memory usage in fast mode
<p align="center">
<img src="https://user-images.githubusercontent.com/71189947/190192041-2315610f-ffa7-4cf9-b6e3-6e2b5b1d3623.png" width="700" height="500"/>
</p>

* At the beggining, the file (or files if `--all`) is loaded into memory as a dataframe : `Loading data in memory` step
* The dataframe is converted into a dictionary, a python object that consums less RAM memory : `Converting dataframe into dictionary` step
* When the conversion is done, the RAM memory occupied by the dataframe is freed : `Memory used by dataframe is freed` step
* Then, the whole computation is done only with the dictionary : `Using dictionary format to work` step

If `--all` is activated, there will be 5 dataframes and 5 dictionaries (not necessarily loaded at the same time, some chromosome files are faster to process like chr2) which use a lot of RAM memory to be 4x faster than slow mode.

### Time comparison
<p align="center">
<img src="https://user-images.githubusercontent.com/71189947/190188112-dae6ee8a-211f-4118-a64f-2f35eb4d692e.png" width="700" height="500"/>
</p>
Since chr1 is the largest file, it limits the computation time. So if the chr1 file takes 2500 sec to be done, then the whole genome should take also around 2500 sec in slow mode for Col-0 (50X). Therefore, the whole genome should take around 600sec in fast mode (be careful with the RAM memory usage).

# Known errors

## IndexError: list index out of range
```bash

Error :

joblib.externals.loky.process_executor._RemoteTraceback:
"""
Traceback (most recent call last):
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/joblib/externals/loky/process_executor.py", line 436, in _process_worker
    r = call_item()
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/joblib/externals/loky/process_executor.py", line 288, in __call__
    return self.fn(*self.args, **self.kwargs)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/joblib/_parallel_backends.py", line 595, in __call__
    return self.func(*args, **kwargs)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/joblib/parallel.py", line 262, in __call__
    return [func(*args, **kwargs)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/joblib/parallel.py", line 262, in <listcomp>
    return [func(*args, **kwargs)
  File "/mnt/data2/rradjas/scripts/ont_to_bam.py", line 267, in main_slow
    create_sam_slow(bedfile_path, SAM_output, BAM_output, genome, args, chrom_name)
  File "/mnt/data2/rradjas/scripts/ont_to_bam.py", line 182, in create_sam_slow
    seq = "".join(change_base(seq, pos_to_change, base))
  File "/mnt/data2/rradjas/scripts/ont_to_bam.py", line 120, in change_base
    if sequence[i] in ["T","A"]:
IndexError: list index out of range
```

Solution : 
Sort the file to avoid having position reversed in `- strand` with :
```bash
sort -k1,1 -k4,4 -k2,2n chr1_methylation.bed > chr1_methylation.sort.bed
```
or to sort all chromosome files :
```bash
for i in {1..5}; do
  sort -k4,4 -k2,2n chr${i}/chr${i}_methylation.bed > chr${i}/chr${i}_methylation.sort.bed
  #rm chr${i}/chr${i}_methylation.bed
done 
```
You can also remove the old file if everything went well.

## KeyError: 'Chr5'

```bash
Running slow mode
Reading genome...
Done
Traceback (most recent call last):
  File "/mnt/data5/rradjas/QLAB/scripts/ont_to_bam.py", line 280, in <module>
    run_function(bedfile, SAM_output, BAM_output, args.genome, args)
  File "/mnt/data5/rradjas/QLAB/scripts/ont_to_bam.py", line 248, in main_slow
    create_sam_slow(bedfile_path, SAM_output, BAM_output, genome, args, chrom_name)
  File "/mnt/data5/rradjas/QLAB/scripts/ont_to_bam.py", line 162, in create_sam_slow
    seq = genome[prev_chrom][first_position:last_position]
KeyError: 'Chr5'
```
You need to use `--all_in_one` to specify that in your file, there is all chromosomes.

## ValueError: invalid literal for int() with base 10: 'C'

```bash
Running fast mode
Reading genome...
Done
 :: Loading dataframe...
 :: Dataframe loaded in 1.46s
 :: Start memory optimisation
Traceback (most recent call last):
  File "/mnt/data5/rradjas/QLAB/scripts/ont_to_bam.py", line 280, in <module>
    run_function(bedfile, SAM_output, BAM_output, args.genome, args)
  File "/mnt/data5/rradjas/QLAB/scripts/ont_to_bam.py", line 235, in main_fast
    df_dict = get_dict(bedfile_path, chrom_name)
  File "/mnt/data5/rradjas/QLAB/scripts/ont_to_bam.py", line 89, in get_dict
    meth_dtf["chr"] = meth_dtf["chr"].apply(lambda chrom: chrom[-1]).astype(np.int8)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/pandas/core/generic.py", line 5912, in astype
    new_data = self._mgr.astype(dtype=dtype, copy=copy, errors=errors)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/pandas/core/internals/managers.py", line 419, in astype
    return self.apply("astype", dtype=dtype, copy=copy, errors=errors)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/pandas/core/internals/managers.py", line 304, in apply
    applied = getattr(b, f)(**kwargs)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/pandas/core/internals/blocks.py", line 580, in astype
    new_values = astype_array_safe(values, dtype, copy=copy, errors=errors)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/pandas/core/dtypes/cast.py", line 1292, in astype_array_safe
    new_values = astype_array(values, dtype, copy=copy)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/pandas/core/dtypes/cast.py", line 1237, in astype_array
    values = astype_nansafe(values, dtype, copy=copy)
  File "/mnt/data2/rradjas/.conda/envs/TE/lib/python3.10/site-packages/pandas/core/dtypes/cast.py", line 1154, in astype_nansafe
    return lib.astype_intsafe(arr, dtype)
  File "pandas/_libs/lib.pyx", line 668, in pandas._libs
  ValueError: invalid literal for int() with base 10: 'C'
 ```
 
You need to remove lines from ChrC and ChrM.

## Packages
Polars : 
```bash
ModuleNotFoundError: No module named 'polars'
```

```bash
Solution :  pip install polars
```
TypeGuard :
```bash
ImportError : cannot import name 'TypeGuard' 
```
```bash
Solution : pip install typing-extensions --upgrade
```

## Custom chromosome name

Try to use conventional names of chromosomes : Chr1 Chr2 Chr3 Chr4 Chr5; in fasta file.
Do not use for example `Chr1_RagTag` name which comes from a genome assembly
