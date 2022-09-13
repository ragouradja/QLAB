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



##############################################################
############              CYTOSINES               ############
##############################################################


#######################################
#######       PROCESSING        #######
#######################################

# Converting from raw multi format read to single format read
multi_to_single_fast5 -i $raw_mutli_fast5_C -s $single_fast5 -t 60 --recursive 

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


# Tombo preprocess
tombo preprocess annotate_raw_with_fastqs --fast5-basedir  $single_fast5 \
--fastq-filenames megalodon_output_C/basecalls.fastq    \
--basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template \
--processes 40 --overwrite 

# Tombo resquiggle
tombo resquiggle $single_fast5 $ref_genome  \
--processes 40 --corrected-group RawGenomeCorrected_000 \
--basecall-group Basecall_1D_000 --overwrite


#######################################
#######       EXTRACTION        #######
#######################################

mkdir datasets -p
mkdir datasets/CG datasets/CHG datasets/CHH -p

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

# Selecting 10M lines
shuf -n 10000000 datasets/CG/samples_CG_poses_positive.tsv > datasets/CG/samples_CG_poses_positive.10m.tsv
# Tar using multiprocess (faster)
pigz datasets/CG/samples_CG_poses_positive.tsv



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

# CHG
# POSITIVE          
deepsignal_plant extract --fast5_dir $single_fast5 \
                         --reference_path $ref_genome \
                         --write_path datasets/CHG/samples_CHG_poses_positive.tsv \
                         --positions $high_conf_positive_CHG \
                         --methy_label 1 \
                         --motifs CHG \
                         --mod_loc 0 \
                         --nproc 20 

shuf -n 10000000 datasets/CHG/samples_CHG_poses_positive.tsv > datasets/CHG/samples_CHG_poses_positive.10m.tsv
pigz datasets/CHG/samples_CG_poses_positive.tsv

# NEGATIVE
deepsignal_plant extract --fast5_dir $single_fast5 \
                         --reference_path $ref_genome \
                         --write_path datasets/CHG/samples_CHG_poses_negative.tsv \
                         --positions $high_conf_negative_CHG \
                         --methy_label 0 \
                         --motifs CHG \
                         --mod_loc 0 \
                         --nproc 20 

shuf -n 10000000 datasets/CHG/samples_CHG_poses_negative.tsv > datasets/CHG/samples_CHG_poses_negative.10m.tsv
pigz datasets/CHG/samples_CHG_poses_negative.tsv


# CHH
# POSITIVE          
deepsignal_plant extract --fast5_dir $single_fast5 \
                         --reference_path $ref_genome \
                         --write_path datasets/CHH/samples_CHH_poses_positive.tsv \
                         --positions $high_conf_positive_CHH \
                         --methy_label 1 \
                         --motifs CHH \
                         --mod_loc 0 \
                         --nproc 20 

shuf -n 10000000 datasets/CHH/samples_CHH_poses_positive.tsv > datasets/CHH/samples_CHH_poses_positive.10m.tsv
pigz datasets/CHH/samples_CHH_poses_positive.tsv

# NEGATIVE
deepsignal_plant extract --fast5_dir $single_fast5 \
                         --reference_path $ref_genome \
                         --write_path datasets/CHH/samples_CHH_poses_negative.tsv \
                         --positions $high_conf_negative_CHH \
                         --methy_label 0 \
                         --motifs CHH \
                         --mod_loc 0 \
                         --nproc 20 

shuf -n 10000000 datasets/CHH/samples_CHH_poses_negative.tsv > datasets/CHH/samples_CHH_poses_negative.10m.tsv
pigz datasets/CHH/samples_CHH_poses_negative.tsv


#######################################
#######       BALANCING         #######
#######################################

# Will create datasets/CHG/balanced/samples_CHG_poses_negative.10m.balanced.tsv
python ${script_dir}/balance_pos_neg.py --pos_file datasets/CHG/samples_CHG_poses_positive.10m.tsv --neg_file datasets/CHG/samples_CHG_poses_negative.10m.tsv --context CHG
pigz datasets/CHG/samples_CHG_poses_negative.10m.tsv

# Will create datasets/CHH/balanced/samples_CHH_poses_negative.10m.balanced.tsv
python ${script_dir}/balance_pos_neg.py --pos_file datasets/CHH/samples_CHH_poses_positive.10m.tsv --neg_file datasets/CHH/samples_CHH_poses_negative.10m.tsv --context CHH
pigz datasets/CHH/samples_CHH_poses_negative.10m.tsv






##############################################################
############                 m6A                  ############
##############################################################

mkdir datasets/A -p
#######################################
#######       PROCESSING        #######
#######################################


##############################
####      POSITIVE        ####
##############################

# Converting from raw multi format read to single format read
multi_to_single_fast5 -i $m6a_raw_mutli_fast5 -s $m6a_single_pos -t 60 --recursive 

# Basecalling using megalodon

megalodon  $m6a_single_pos  \
--guppy-params "-d $rerio_model --num_callers 8 --gpu_runners_per_device 64 --chunks_per_runner 256 --chunk_size 2000" \
--guppy-server-path /usr/bin/guppy_basecall_server \
--guppy-config  res_dna_r941_min_modbases-all-context_v001.cfg  \
--guppy-timeout 120 \
--outputs basecalls \
--output-directory megalodon_output_A \
--reference $ref_genome \
--devices 0 --processes 30 --overwrite


# Tombo preprocess
tombo preprocess annotate_raw_with_fastqs --fast5-basedir  $m6a_single_pos \
--fastq-filenames megalodon_output_A/basecalls.fastq    \
--basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template \
--processes 40 --overwrite 

# Tombo resquiggle
tombo resquiggle $m6a_single_pos $ref_genome  \
--processes 40 --corrected-group RawGenomeCorrected_000 \
--basecall-group Basecall_1D_000 --overwrite




#######################################
#######       EXTRACTION        #######
#######################################



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




##############################################################
############          MODEL DATASETS              ############
##############################################################

mkdir training -p
mkdir training/train training/valid training/model -p

#######################################
#######   TRAIN/VALIDATION SET  #######
#######################################

# CONCATENATION OF ALL DATA INTO ONE FILE
cat datasets/A/samples_A_positive.10m.tsv  \
    datasets/A/samples_A_negative.10m.tsv \
    datasets/CG/samples_CG_poses_positive.10m.tsv \
    datasets/CG/samples_CG_poses_negative.10m.tsv \
    datasets/CHG/samples_CHG_poses_positive.10m.tsv \
    datasets/CHG/balanced/samples_CHG_poses_negative.10m.balanced.tsv \
    datasets/CHH/samples_CHH_poses_positive.10m.tsv \
    datasets/CHH/balanced/samples_CHH_poses_negative.10m.balanced.tsv | shuf > training/concat_A_CG_CHG_CHH.80m.tsv

# TRAIN SET (90% of data)
head -n 72000000  training/concat_A_CG_CHG_CHH.80m.tsv > training/train/samples_A_CG_CHG_CHH.72m.train.tsv
# VALIDATION SET (10% of data)
tail -n 8000000  training/concat_A_CG_CHG_CHH.80m.tsv > training/valid/samples_A_CG_CHG_CHH.8m.valid.tsv

rm training/concat_A_CG_CHG_CHH.80m.tsv



##############################################################
############              TRAINING                ############
##############################################################

CUDA_VISIBLE_DEVICES=0 deepsignal_plant train --train_file training/train/samples_A_CG_CHG_CHH.72m.train.tsv \
                                              --valid_file training/valid/samples_A_CG_CHG_CHH.8m.valid.tsv \
                                              --model_dir training/model/model_10M \
                                            #   --batch_size 2048


