# Stretch detection from DSP tsv file
Decription of stretch
image ?

To be launch after DSP pipeline and the tsv_to_bed script.

## How it works
This bash script will first create some intermediate files and folders to optimize the calculations. 
The steps are : 
* Extract each chromosoms and sort the lines by : Chr, read name, start position
The sort is needed to avoid having decreasing start position in Crick strand.
```bash
sort -T tmp -k1,1 -k4,4 -k2,2n
```


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
│   ├── chr2_methylation.bed
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
│   ├── chr3_methylation.bed
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
│   ├── chr4_methylation.bed
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
│   ├── chr5_methylation.bed
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

## How to run
## Outputs
## Known errors
