tsv_file=/mnt/data2/rradjas/ONT/Rdr2/dsp_rdr2/tsv/fast5s.C.call_mods.tsv.gz
reads_analysis_folder=/mnt/data2/rradjas/ONT/Rdr2/dsp_rdr2/tsv/reads_analysis
scripts_folder=/mnt/data2/rradjas/scripts
genome_ref=/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta


get_stretch(){
    local CHR=$1
    chr=${CHR,,} # Chr1 --> chr1
    mkdir ${reads_analysis_folder}/${chr}/stretch/pvalue/ -p
    mkdir ${reads_analysis_folder}/${chr}/context/freq -p
    echo "Doing " $CHR
    for context in CG CHG CHH; do
        echo "Getting methylation frequencies for " ${context} context
        python ${scripts_folder}/correct_freq.py --file $reads_analysis_folder/${chr}/context/${context}_context.bed \
                                                 --output_file $reads_analysis_folder/${chr}/context/freq/${context}_freq.bed

        echo "Sorting" ${context} context
        sort -T tmp -k1,1 -k4,4 -k2,2n $reads_analysis_folder/${chr}/context/${context}_context.bed > $reads_analysis_folder/${chr}/context/${context}_context.sort.bed
        rm  $reads_analysis_folder/${chr}/context/${context}_context.bed 

        echo "Compute stretches detection" ${context} context
        python ${scripts_folder}/get_stretch.py --read_file $reads_analysis_folder/${chr}/context/${context}_context.sort.bed \
                                                --meth_file $reads_analysis_folder/${chr}/context/freq/${context}_freq.bed \
                                                --output_file $reads_analysis_folder/${chr}/stretch/pvalue/${context}_stretch_pvalue.bed

        echo "Done for " ${context} context
    done
    echo "Done for " $CHR
}


echo "Getting chr and context files"
python ${scripts_folder}/separator.py --file $tsv_file \
                                      -o $reads_analysis_folder \
                                      --tsv  \
                                      -s 3
mkdir tmp -p
start_time=$(date +%s) 


for CHR in Chr1 Chr2 Chr3 Chr4 Chr5; do

    # DOING ALL CHROMOSOMS IN PARALLEL : relatively high memory usage
    get_stretch $CHR & 

    # DOING ONE BY ONE : less memory usage
    # get_stretch $CHR   
done
wait


mkdir ${reads_analysis_folder}/genome_stretch/ -p
mkdir ${reads_analysis_folder}/bed/ -p
for context in CG CHG CHH; do
cat ${reads_analysis_folder}/chr*/stretch/pvalue/${context}_stretch_pvalue.bed > ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch.bed
thresh=$(echo "-l(0.05/$(cat ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch.bed | wc -l)) / l(10)" | bc -l)
awk -v Thresh="$thresh" ' $7 >= Thresh' ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch.bed > ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch_sign.bed

# cat ${reads_analysis_folder}/chr*/context/freq/${context}_freq.bed > ${reads_analysis_folder}/bed/${context}_freq.bed 

done;



end_time=$(date +%s) 

echo "Total time : " $(( end_time - start_time ))
