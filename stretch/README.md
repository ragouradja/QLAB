# Stretch detection from DSP tsv file
Decription of stretch
image ?

To be launch after DSP pipeline.

## How it works
This bash script will first create some intermediate files and folders to optimize the calculations.


<details>
  
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
```
  
</details>

```bash
sort -T tmp -k1,1 -k4,4 -k2,2n
```
## How to run
## Outputs
## Known errors
