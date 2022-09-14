# Train a deep learning model using DSP



# How it works

## 1. Path variables


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

* raw_mutli_fast5_C

Path to the folder containing raw fast5 with cytosines sequencing data in multi-read format. These fast5 are going to be converted in single-read format.

If the fast5 are already in single-read format, this variable can be empty (in this case, multi_to_single_fast5 will throw an error that can be ignored)

* single_fast5 : Path to the folder that will contain fast5 in single-read format.

* high_conf_positive_C* : File with positions of high confidence positive cytosines having methylation level >= 0.9 and >= 5X from BS or / and ONT.

* high_conf_negative_C* : File with positions of high confidence negative cytosines having methylation level == 0 and >= 5X from BS or / and ONT.

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


* m6a_raw_mutli_fast5 : Same as raw_mutli_fast5_C but for positive m6A.

* m6a_single_pos : Same as single_fast5 but for positive m6A. Negative data for m6A are extracted from single_fast5.

* ref_genome : Path to the reference genome.

* script_dir : Folder with balance_pos_neg.py


## 2. Cytosines data
### 2.1 Processing steps

The raw fast5 reads are converted from multi-read format to single-read format. If your fast5 are already in multi-read format, *$raw_mutli_fast5_C* can be empty to skip this step.

```bash
multi_to_single_fast5 -i $raw_mutli_fast5_C -s $single_fast5 -t 60 --recursive 
```

Then, basecalling the reads using Megalodon to have the fastq file.

```bash
# Basecalling using megalodon

megalodon  $single_fast5  \
--guppy-params "-d $rerio_model --num_callers 8 --gpu_runners_per_device 64 --chunks_per_runner 256 --chunk_size 2000" \
--guppy-server-path /usr/bin/guppy_basecall_server \
--guppy-config  res_dna_r941_min_modbases_5mC_v001.cfg  \
--guppy-timeout 120 \
--outputs basecalls \
--output-directory megalodon_output_C \
--reference $ref_genome \
--devices 0 --processes 30 --overwrite
```

Finally, preprocessing and resquiggle step with Tombo.

```bash
# Tombo preprocess
tombo preprocess annotate_raw_with_fastqs --fast5-basedir  $single_fast5 \
--fastq-filenames megalodon_output_C/basecalls.fastq    \
--basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template \
--processes 40 --overwrite 

# Tombo resquiggle
tombo resquiggle $single_fast5 $ref_genome  \
--processes 40 --corrected-group RawGenomeCorrected_000 \
--basecall-group Basecall_1D_000 --overwrite

```

### 2.2 Data extraction

Extraction example for positive CG data using high confidence positions of methylated cytosines :
```bash

# CG
# POSITIVE          
deepsignal_plant extract --fast5_dir $single_fast5 \
                         --reference_path $ref_genome \
                         --write_path datasets/CG/samples_CG_poses_positive.tsv \
                         --positions $high_conf_positive_CG \
                         --methy_label 1 \
                         --motifs CG \
                         --mod_loc 0 \
                         --nproc 20 
```

After extraction, you will have a certain amount of lines (cytosines), **make sure that you have at least 10M lines if you want to create datasets with 10M positive and 10M negative.**

```bash
# Randomly selecting 10M lines
shuf -n 10000000 datasets/CG/samples_CG_poses_positive.tsv > datasets/CG/samples_CG_poses_positive.10m.tsv
```
Then, using pigz to compress the extraction file.

```bash
# Tar using multiprocess (faster)
pigz datasets/CG/samples_CG_poses_positive.tsv
```

Extraction of negative CG data using high confidence positions of unmethylated cytosines :
```bash
# NEGATIVE
deepsignal_plant extract --fast5_dir $single_fast5 \
                         --reference_path $ref_genome \
                         --write_path datasets/CG/samples_CG_poses_negative.tsv \
                         --positions $high_conf_negative_CG \
                         --methy_label 0 \
                         --motifs CG \
                         --mod_loc 0 \
                         --nproc 20 

shuf -n 10000000 datasets/CG/samples_CG_poses_negative.tsv > datasets/CG/samples_CG_poses_negative.10m.tsv
pigz datasets/CG/samples_CG_poses_negative.tsv
```
For this extraction, you will have for sure more than 10M lines (much more) since there is more unmethylated cytosines than methylated.

### 2.3 Balancing step for CHG and CHH
Balancing kmers of methylated and unmethylated cytosines to have similar proportions of kmers (see DSP paper). This will create a new negative dataset to use. The old one is compressed.

```bash
# Will create datasets/CHG/balanced/samples_CHG_poses_negative.10m.balanced.tsv
python ${script_dir}/balance_pos_neg.py --pos_file datasets/CHG/samples_CHG_poses_positive.10m.tsv 
                                        --neg_file datasets/CHG/samples_CHG_poses_negative.10m.tsv 
                                        --context CHG
pigz datasets/CHG/samples_CHG_poses_negative.10m.tsv

# Will create datasets/CHH/balanced/samples_CHH_poses_negative.10m.balanced.tsv
python ${script_dir}/balance_pos_neg.py --pos_file datasets/CHH/samples_CHH_poses_positive.10m.tsv
                                        --neg_file datasets/CHH/samples_CHH_poses_negative.10m.tsv 
                                        --context CHH
pigz datasets/CHH/samples_CHH_poses_negative.10m.tsv


```

## 3. m6A data
### 3.1 Processing steps
Same as 2.1.

### 3.2 Data extraction

For positive m6a, there is no file with high confidence position (or you can add one). 

```bash
# A
# POSITIVE          
deepsignal_plant extract --fast5_dir $m6a_single_pos \
                         --reference_path $ref_genome \
                         --write_path datasets/A/samples_A_positive.tsv \
                         --methy_label 1 \
                         --motifs A \
                         --mod_loc 0 \
                         --nproc 20 

# Selecting 10M lines
shuf -n 10000000 datasets/A/samples_A_positive.tsv  > datasets/A/samples_A_positive.10m.tsv 
# Tar using multiprocess (faster)
pigz datasets/A/samples_A_positive.tsv 
```

The negative m6a data are extracted from fast5 reads used to extract cytosines data. If you don't want to use these fast5, you can change *$single_fast5* by a new variable (*$m6a_single_neg* for example) containing path to folder with fast5 to use to extract negative m6a. Make sure that these fast5 have been preprocessed (2.1 Processing steps)

```bash
# NEGATIVE          
deepsignal_plant extract --fast5_dir $single_fast5 \
                         --reference_path $ref_genome \
                         --write_path datasets/A/samples_A_negative.tsv \
                         --methy_label 0 \
                         --motifs A \
                         --mod_loc 0 \
                         --nproc 20 

# Selecting 10M lines
shuf -n 10000000 datasets/A/samples_A_negative.tsv  > datasets/A/samples_A_negative.10m.tsv 
# Tar using multiprocess (faster)
pigz datasets/A/samples_A_negative.tsv 

```


## 4. Training datasets


### 4.1 Concatenation of all datasets
Now that we have a dataset for each context of 10M positive + 10M negative data, we can concatenate (and shuffle) them all to create a file with 80M. 

```bash
cat datasets/A/samples_A_positive.10m.tsv  \
    datasets/A/samples_A_negative.10m.tsv \
    datasets/CG/samples_CG_poses_positive.10m.tsv \
    datasets/CG/samples_CG_poses_negative.10m.tsv \
    datasets/CHG/samples_CHG_poses_positive.10m.tsv \
    datasets/CHG/balanced/samples_CHG_poses_negative.10m.balanced.tsv \
    datasets/CHH/samples_CHH_poses_positive.10m.tsv \
    datasets/CHH/balanced/samples_CHH_poses_negative.10m.balanced.tsv | shuf > training/concat_A_CG_CHG_CHH.80m.tsv
   
```

### 4.1 Training set
From this file, we create the training set and validation set by selecting 90% and 10% of these 80M (you can change these proportion). These sets are going to be used directly in the training process.

```bash
# TRAIN SET (90% of data)
head -n 72000000  training/concat_A_CG_CHG_CHH.80m.tsv > training/train/samples_A_CG_CHG_CHH.72m.train.tsv

```
### 4.2 Validation set
```bash
# VALIDATION SET (10% of data)
tail -n 8000000  training/concat_A_CG_CHG_CHH.80m.tsv > training/valid/samples_A_CG_CHG_CHH.8m.valid.tsv

rm training/concat_A_CG_CHG_CHH.80m.tsv

```

## 5. Model training

```bash
CUDA_VISIBLE_DEVICES=0 deepsignal_plant train --train_file training/train/samples_A_CG_CHG_CHH.72m.train.tsv \
                                              --valid_file training/valid/samples_A_CG_CHG_CHH.8m.valid.tsv \
                                              --model_dir training/model/model_10M \
                                            #   --batch_size 2048
```

