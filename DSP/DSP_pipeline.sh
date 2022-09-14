fast5_dir=/mnt/data2/rradjas/ONT/Rdr2/fast5_single
genome_ref=/mnt/data2/rradjas/assembly/Col-CEN/Col-CEN_v1.2.fasta 
megalodon_output=megalodon_output_rdr2
fastq_file=/mnt/data2/rradjas/ONT/Rdr2/${megalodon_output}/basecalls.fastq
output_dir=/mnt/data2/rradjas/ONT/Rdr2/dsp_rdr2
raw_fast5=/mnt/data2/rradjas/ONT/Rdr2/raw

multi_to_single_fast5 -i $raw_fast5 -s $fast5_dir -t 90 --recursive
echo END_MULTI

rerio_model=/users/a2e/quadrana/miniconda3/lib/python3.8/site-packages/megalodon/rerio/basecall_models/
megalodon  $fast5_dir \
--guppy-params "-d $rerio_model --num_callers 8 --gpu_runners_per_device 64 --chunks_per_runner 256 --chunk_size 2000" \
--guppy-server-path /usr/bin/guppy_basecall_server \
--guppy-config  res_dna_r941_min_modbases_5mC_v001.cfg  \
--guppy-timeout 120 \
--outputs basecalls \
--output-directory ${megalodon_output} \
--reference /mnt/data2/rradjas/assembly/Col-CEN/Col-CEN_v1.2.fasta \
--devices 0 --processes 30 --overwrite
echo END_BASECALLING


tombo preprocess annotate_raw_with_fastqs --fast5-basedir  $fast5_dir \
--fastq-filenames $fastq_file \
--basecall-group Basecall_1D_000 --basecall-subgroup BaseCalled_template \
--processes 20 --overwrite 
echo END_TOMBOPRE


tombo resquiggle $fast5_dir  $genome_ref \
--processes 20 --corrected-group RawGenomeCorrected_000 \
--basecall-group Basecall_1D_000 --overwrite 
echo END_TOMBORES


mkdir $output_dir -p
deepsignal_plant call_mods --input_path $fast5_dir \
--model_path model.dp2.CNN.arabnrice2-1_120m_R9.4plus_tem.bn13_sn16.both_bilstm.epoch6.ckpt \
--result_file $output_dir/fast5s.C.call_mods.tsv \
--corrected_group RawGenomeCorrected_000 \
--reference_path $genome_ref \
--motifs C --nproc 20 --nproc_gpu 5
echo END_DSPMODS
