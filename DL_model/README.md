# Train a deep learning model using DSP



# How it work

## 1. Path variables
* raw_mutli_fast5_C

Path to the folder containing raw fast5 with cytosines sequencing data in multi-read format. These fast5 are going to be converted in single-read format.

If the fast5 are already in single-read format, this variable can be empty (in this case, multi_to_single_fast5 will throw an error that can be ignored)

* single_fast5
Path to the folder that will contain fast5 in single-read format.

* high_conf_positive_C*
File with positions of high confidence positive cytosines having methylation level >= 0.9 and >= 5X from BS or / and ONT.

* high_conf_negative_C*
File with positions of high confidence negative cytosines having methylation level == 0 and >= 5X from BS or / and ONT.

Example of position file with columns Chr end strand: 

```bash
Chr1    3724    -
Chr1    3804    -
Chr1    3958    -
Chr1    3963    +
Chr1    3964    -
Chr1    4490    +
Chr1    4639    -
Chr1    4658    -
Chr1    4688    -
Chr1    9323    -
```


```bash
raw_mutli_fast5_C=/mnt/data5/rradjas/all_contexts/m6a/test_script/neg
single_fast5=single_C

high_conf_positive_CG=/mnt/data5/rradjas/all_contexts/BS_data/CG/poses/methylated_09_5_biseq_CG.bed
high_conf_negative_CG=/mnt/data5/rradjas/all_contexts/BS_data/CG/poses/unmethylated_0_5_biseq_CG.bed

high_conf_positive_CHG=/mnt/data5/rradjas/all_contexts/BS_data/CHG/poses/methylated_09_5_biseq_CHG.bed
high_conf_negative_CHG=/mnt/data5/rradjas/all_contexts/BS_data/CHG/poses/unmethylated_0_5_biseq_CHG.bed

high_conf_positive_CHH=/mnt/data5/rradjas/all_contexts/BS_data/CHH/poses/methylated_09_5_biseq_CHH.bed
high_conf_negative_CHH=/mnt/data5/rradjas/all_contexts/BS_data/CHH/poses/unmethylated_0_5_biseq_CHH.bed


m6a_raw_mutli_fast5=/mnt/data5/rradjas/all_contexts/m6a/test_script/pos
m6a_single_pos=single_A
ref_genome=/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_split.fasta

script_dir=/mnt/data2/rradjas/scripts
rerio_model=/users/a2e/quadrana/miniconda3/lib/python3.8/site-packages/megalodon/rerio/basecall_models/
```

