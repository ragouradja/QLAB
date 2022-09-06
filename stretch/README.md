# Stretch detection from DSP tsv file
Decription of stretch
image ?

To be launch after DSP pipeline and the tsv_to_bed script.


```bash
Chr1    16501609        16501610        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501603        16501604        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501600        16501601        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501598        16501599        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501597        16501598        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501592        16501593        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501590        16501591        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501586        16501587        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501585        16501586        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
Chr1    16501582        16501583        073be7c9-69c5-498d-b8c7-b28a6124ca5c    -       0
```



## How to run


Make sure that you put the correct paths for these 4 variables before running the script.
```bash
bed_file=/mnt/data2/rradjas/ONT/Rdr2/dsp_rdr2/tsv/methylation.bed
reads_analysis_folder=/mnt/data2/rradjas/ONT/Rdr2/dsp_rdr2/tsv/reads_analysis

scripts_folder=/mnt/data2/rradjas/scripts
genome_ref=/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta
```

```bash
bash stretch_pipeline.sh
```

## How it works
This bash script will first create some intermediate files and folders to optimize the calculations. 
The steps are : 
### 1. Extract and sort
Each chromosom are extracted from the inital bedfile and sorted by : read name, start position. The sort is needed to avoid having decreasing start position in - strand.
```bash
grep $CHR $bed_file | sort -T tmp -k4,4 -k2,2n > ${reads_analysis_folder}/${chr}/${chr}_methylation.sort.bed
```
#### Output : chrX/chrX_methylation.sort.bed

### 2. Context separation
Each chromosom file is separated into CG, CHG and CHH context files. 
```bash
python ${scripts_folder}/sep_context.py  --file ${reads_analysis_folder}/${chr}/${chr}_methylation.sort.bed \
                                         --genome $genome_ref \
                                         --output_dir $reads_analysis_folder/${chr}/context \
                                         -s 5
```
#### Output : chrX/context/CG_context.bed

### 3. Methylation frequencies
From each context files, methylation frequencies are computed to create a bedfile.
```bash
python ${scripts_folder}/correct_freq.py --file $reads_analysis_folder/${chr}/context/${context}_context.bed \
                                         --output_file $reads_analysis_folder/${chr}/context/freq/${context}_freq.bed
```
At the end, all these freq files from all chromosoms are concatened into one file per context with genome wide information.
#### Output : chrX/context/freq/CG_freq.bed

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


### 4. Concatenating final outputs
These stretch file per chromosom per context are concatenated into one file per context having genome wide stretches information.

```bash
mkdir ${reads_analysis_folder}/genome_stretch/ -p
mkdir ${reads_analysis_folder}/bed/ -p
for context in CG CHG CHH;
  cat ${reads_analysis_folder}/chr*/stretch/pvalue/${context}_stretch_pvalue_sign.bed > ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch_sign.bed
  cat ${reads_analysis_folder}/chr*/context/freq/${context}_freq.bed > ${reads_analysis_folder}/bed/${context}_freq.bed 
done;
```
#### Output : 
* genome_stretch/CG_genome_stretch_sign.bed
* bed/CG_freq.bed



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

## Outputs
## Known errors
