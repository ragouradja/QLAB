tsv_file=/mnt/data5/rradjas/ONT/Col-0/dsp_col-0/tsv/
reads_analysis_folder=/mnt/data5/rradjas/ONT/Col-0/dsp_col-0/reads_analysis

scripts_folder=/mnt/data2/rradjas/scripts
genome_ref=/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta


get_stretch(){
    local CHR=$1
    chr=${CHR,,} # Chr1 --> chr1
    mkdir ${reads_analysis_folder}/${chr}/stretch/pvalue/ -p
    mkdir ${reads_analysis_folder}/${chr}/context/freq -p
    echo "Doing " $CHR
    
    # Comment these 2 lines if you have already your methylation bedfile created for each chr in the reads_analysis folder
    echo "Grep and sort" $CHR
    grep $CHR $read_file | sort -T tmp -k1,1 -k4,4 -k2,2n > ${reads_analysis_folder}/${chr}_methylation.sort.bed

    echo "Context separation"
    python ${scripts_folder}/sep_context.py  --file ${reads_analysis_folder}/${chr}/${chr}_methylation.bed --genome $genome_ref --output_dir $reads_analysis_folder/${chr}/context

    for context in CG CHG CHH; do
        echo "Getting methylation frequencies for " ${context} context
        python ${scripts_folder}/correct_freq.py $reads_analysis_folder/${chr}/context/${context}_context.bed $reads_analysis_folder/${chr}/context/freq/${context}_freq.bed

        # echo "Sorting" ${context} context
        # sort -T tmp -k1,1 -k4,4 -k2,2n $reads_analysis_folder/${chr}/context/${context}_context.bed > $reads_analysis_folder/${chr}/context/${context}_context.sort.bed
        # rm  $reads_analysis_folder/${chr}/context/${context}_context.bed 

        echo "Compute stretches detection" ${context} context
        python ${scripts_folder}/get_stretch.py $reads_analysis_folder/${chr}/context/${context}_context.sort.bed \
                $reads_analysis_folder/${chr}/context/freq/${context}_freq.bed \
                $reads_analysis_folder/${chr}/stretch/pvalue/${context}_stretch_pvalue.bed

        thresh=$(echo "-l(0.05/$(cat ${reads_analysis_folder}/${chr}/stretch/pvalue/${context}_stretch_pvalue.bed  | wc -l)) / l(10)" | bc -l)
        awk -v Thresh="$thresh" ' $7 >= Thresh' ${reads_analysis_folder}/${chr}/stretch/pvalue/${context}_stretch_pvalue.bed  > ${reads_analysis_folder}/${chr}/stretch/pvalue/${context}_stretch_pvalue_sign.bed
        

        echo "Done for " ${context} context
    done
    echo "Done for " $CHR
}




start_time=$(date +%s) 


for CHR in Chr1 Chr2 Chr3 Chr4 Chr5; do

    # DOING ALL CHROMOSOMS IN PARALLEL : relatively high memory usage
    # get_stretch $CHR & 

    # DOING ONE BY ONE : less memory usage
    get_stretch $CHR   
done
wait


mkdir ${reads_analysis_folder}/genome_stretch/ -p
mkdir ${reads_analysis_folder}/bed/ -p
for context in CG CHG CHH;
cat ${reads_analysis_folder}/chr*/stretch/pvalue/${context}_stretch_pvalue_sign.bed > ${reads_analysis_folder}/genome_stretch/${context}_genome_stretch_sign.bed
cat ${reads_analysis_folder}/chr*/context/freq/${context}_freq.bed > ${reads_analysis_folder}/bed/${context}_freq.bed 
done;


end_time=$(date +%s) 

echo "Total time : " $(( end_time - start_time ))