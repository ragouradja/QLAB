# Stretch detection from DSP tsv file


**To be launch after DSP pipeline**

## How to run


Make sure that you put the correct paths for these 4 variables before running the script.
```bash
tsv_file=/mnt/data2/rradjas/ONT/Rdr2/dsp_rdr2/tsv/fast5s.C.call_mods.tsv.gz
reads_analysis_folder=/mnt/data2/rradjas/ONT/Rdr2/dsp_rdr2/tsv/reads_analysis
scripts_folder=/mnt/data2/rradjas/scripts
genome_ref=/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta
```

The tsv_file can be a compressed file ending with "gz".


## How it works
This bash script will first create some intermediate files and folders to optimize the calculations. 
The steps are : 

### 1. Extract and sort
From the inital TSV file, this script will extract all 5 chromosomes and cytosines per context of each chromosome.

`-s` is used to specify the column containing the strand information. `--tsv` indicates that the script will select only required columns for analysis from the TSV file. 

```bash
echo "Getting chr and context files"
python ${scripts_folder}/separator.py --file $tsv_file \
                                      -o $reads_analysis_folder \
                                      --tsv  \
                                      -s 3
```
#### Output : 
* chrX/chrX_methylation.bed
* chrX/context/CG_context.bed

### 2. Methylation frequencies
From each context files, methylation frequencies are computed to create a bedfile.
```bash
python ${scripts_folder}/correct_freq.py --file $reads_analysis_folder/${chr}/context/${context}_context.bed \
                                         --output_file $reads_analysis_folder/${chr}/context/freq/${context}_freq.bed
```
At the end, all these freq files from all chromosoms are concatened into one file per context with genome wide information.
#### Output : chrX/context/freq/CG_freq.bed

### 3. Sorting context files
These files need to be sorted by READ and START_POSITION to avoid having decreasing start position in - strand (which will can be problematic for finding stretch).

```bash
sort -T tmp -k4,4 -k2,2n $reads_analysis_folder/${chr}/context/${context}_context.bed > $reads_analysis_folder/${chr}/context/${context}_context.sort.bed
rm  $reads_analysis_folder/${chr}/context/${context}_context.bed 
```

### 4. Stretch detection
Using the bedfile with read informations to detect consecutively methylated cytosines and using the methylation frequencies per position to compute the pvalue for each stretch.
```bash
python ${scripts_folder}/get_stretch.py --read_file $reads_analysis_folder/${chr}/context/${context}_context.sort.bed \ 
                                        --meth_file $reads_analysis_folder/${chr}/context/freq/${context}_freq.bed \
                                        --output_file $reads_analysis_folder/${chr}/stretch/pvalue/${context}_stretch_pvalue.bed
```
#### Output : chrX/stretch/pvalue/CG_stretch_pvalue.bed

Doing a Bonferroni correction on $alpha$ threshold by dividing $alpha=0.05$ by the total number of stretches in one chromosom for one context and then selecting statistically significant stretches. 
```bash
thresh=$(echo "-l(0.05/$(cat ${reads_analysis_folder}/${chr}/stretch/pvalue/${context}_stretch_pvalue.bed  | wc -l)) / l(10)" | bc -l)
awk -v Thresh="$thresh" ' $7 >= Thresh' ${reads_analysis_folder}/${chr}/stretch/pvalue/${context}_stretch_pvalue.bed  > ${reads_analysis_folder}/${chr}/stretch/pvalue/${context}_stretch_pvalue_sign.bed

```
#### Output : chrX/stretch/pvalue/CG_stretch_pvalue_sign.bed


### 5. Concatenating final outputs
These stretch file per chromosom per context are concatenated into one file per context having genome wide stretches information. Then, statistically significant stretches are extracted using Bonferroni correction.

```bash
mkdir ${reads_analysis_folder}/genome_stretch/ -p
mkdir ${reads_analysis_folder}/bed/ -p

for context in CG CHG CHH; do
cat ${reads_analysis_folder}/chr*/stretch/pvalue/${context}_stretch_pvalue.bed > ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch.bed
thresh=$(echo "-l(0.05/$(cat ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch.bed | wc -l)) / l(10)" | bc -l)
awk -v Thresh="$thresh" ' $7 >= Thresh' ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch.bed > ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch_sign.bed
cat ${reads_analysis_folder}/chr*/context/freq/${context}_freq.bed > ${reads_analysis_folder}/bed/${context}_freq.bed 
done;
```

This command line will concatenate all chromosomes frequencies files into one per context.
```bash
cat ${reads_analysis_folder}/chr*/context/freq/${context}_freq.bed > ${reads_analysis_folder}/bed/${context}_freq.bed 
```


#### Output : 
* genome_stretch/CG_genome_stretch_sign.bed
* bed/CG_freq.bed



### Parallelization

All 5 chromosomes can be processed in parallel using `&`. But if you are short in RAM memory, you can put `#` on the line `get_stretch $CHR & ` and uncomment `get_stretch $CHR ` to process chromosomes files iteratively (will be longer).

```bash
for CHR in Chr1 Chr2 Chr3 Chr4 Chr5; do

    # DOING ALL CHROMOSOMS IN PARALLEL : relatively high memory usage
    get_stretch $CHR & 

    # DOING ONE BY ONE : less memory usage
    # get_stretch $CHR   
done
wait
```
Unlike the `ont_to_bam.py` script where the `fast mode` is loading the whole chromosome in RAM memory, here it is only the context files of each chromosome that will be loaded. The parallelization of these file should use less RAM memory than `ont_to_bam.py` in `fast mode`.

<details>
  <summary>Files and folders created</summary>
  
```bash

reads_analysis/

├── chr1
│   ├── chr1_methylation.sort.bed
│   ├── context
│   │   ├── CG_context.bed
│   │   ├── CHG_context.bed
│   │   ├── CHH_context.bed
│   │   └── freq
│   │       ├── CG_freq.bed
│   │       ├── CHG_freq.bed
│   │       └── CHH_freq.bed
│   ├── stretch
│   │   └── pvalue
│   │       ├── CG_stretch_pvalue.bed
│   │       ├── CG_stretch_pvalue_sign.bed
│   │       ├── CHG_stretch_pvalue.bed
│   │       ├── CHG_stretch_pvalue_sign.bed
│   │       ├── CHH_stretch_pvalue.bed
│   │       └── CHH_stretch_pvalue_sign.bed
├── chr2
│   ├── chr2_methylation.sort.bed
│   ├── context
│   │   ├── CG_context.bed
│   │   ├── CHG_context.bed
│   │   ├── CHH_context.bed
│   │   └── freq
│   │       ├── CG_freq.bed
│   │       ├── CHG_freq.bed
│   │       └── CHH_freq.bed
│   └── stretch
│       └── pvalue
│           ├── CG_stretch_pvalue.bed
│           ├── CG_stretch_pvalue_sign.bed
│           ├── CHG_stretch_pvalue.bed
│           ├── CHG_stretch_pvalue_sign.bed
│           ├── CHH_stretch_pvalue.bed
│           └── CHH_stretch_pvalue_sign.bed
├── chr3
│   ├── chr3_methylation.sort.bed
│   ├── context
│   │   ├── CG_context.bed
│   │   ├── CHG_context.bed
│   │   ├── CHH_context.bed
│   │   └── freq
│   │       ├── CG_freq.bed
│   │       ├── CHG_freq.bed
│   │       └── CHH_freq.bed
│   └── stretch
│       └── pvalue
│           ├── CG_stretch_pvalue.bed
│           ├── CG_stretch_pvalue_sign.bed
│           ├── CHG_stretch_pvalue.bed
│           ├── CHG_stretch_pvalue_sign.bed
│           ├── CHH_stretch_pvalue.bed
│           └── CHH_stretch_pvalue_sign.bed
├── chr4
│   ├── chr4_methylation.sort.bed
│   ├── context
│   │   ├── CG_context.bed
│   │   ├── CHG_context.bed
│   │   ├── CHH_context.bed
│   │   └── freq
│   │       ├── CG_freq.bed
│   │       ├── CHG_freq.bed
│   │       └── CHH_freq.bed
│   └── stretch
│       └── pvalue
│           ├── CG_stretch_pvalue.bed
│           ├── CG_stretch_pvalue_sign.bed
│           ├── CHG_stretch_pvalue.bed
│           ├── CHG_stretch_pvalue_sign.bed
│           ├── CHH_stretch_pvalue.bed
│           └── CHH_stretch_pvalue_sign.bed
├── chr5
│   ├── chr5_methylation.sort.bed
│   ├── context
│   │   ├── CG_context.bed
│   │   ├── CHG_context.bed
│   │   ├── CHH_context.bed
│   │   └── freq
│   │       ├── CG_freq.bed
│   │       ├── CHG_freq.bed
│   │       └── CHH_freq.bed
│   └── stretch
│       └── pvalue
│           ├── CG_stretch_pvalue.bed
│           ├── CG_stretch_pvalue_sign.bed
│           ├── CHG_stretch_pvalue.bed
│           ├── CHG_stretch_pvalue_sign.bed
│           ├── CHH_stretch_pvalue.bed
│           └── CHH_stretch_pvalue_sign.bed
├── genome_stretch
│   ├── CG_genome_stretch_sign.bed
│   ├── CHG_genome_stretch_sign.bed
│   └── CHH_genome_stretch_sign.bed
├── bed
│   ├── CG_freq.bed
│   ├── CHG_freq.bed
│   └── CHH_freq.bed
```
  
</details>
