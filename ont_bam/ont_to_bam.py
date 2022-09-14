""" Create bed file with boundings of reads"""

"""
Sort the file to avoid having position reversed in minus strand with : sort -k1,1 -k4,4 -k2,2n chr1_methylation.bed > chr1_methylation.sort.bed
Error :

###
Chr1 / Chr1_RagTag --> custom name of chr
Custom chr name are not processed correctly for now
###


###
input chrom for all; folder architecture :

../data/chr1/meth_chr1.bed
../data/chr2/meth_chr2.bed
../data/chr3/meth_chr3.bed
../data/chr4/meth_chr4.bed
../data/chr5/meth_chr5.bed

or

../data/chrom_folder/meth_chr1.bed
../data/chrom_folder/meth_chr2.bed
../data/chrom_folder/meth_chr3.bed
../data/chrom_folder/meth_chr4.bed
../data/chrom_folder/meth_chr5.bed
###


###
Polars : 
ImportError : cannot import name 'TypeGuard' 
Solution : pip install typing-extensions --upgrade
###


###
KeyError: '0001e11e-a6bd-4a0b-83b3-989a7bdfe519'
--> Read name not found in fasta file : fastaFromBed didn't work
Solution : Give a genome fasta file with formated lines 

wrong fasta file : 
>Chr1
sequence chr1 all in one line
>Chr2
sequence chr2 all in one line

correct fasta file:
>Chr1
TAAAACCTAAAACCTAAAACCTAAAACCTAAAACCTAAACCCTAAACCCTAAAACCCTAA
ACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAACCTAAACCT 
...
###


"""


"""
BENCHMARK:

27GO bedfile of 50X data : 
Fast ont_to_bam but with RAM memory usage: 1258s
Slow ont_to_bam but few RAM memory usage : 3000s
FASTER : 432s

"""


import pandas as pd
import numpy as np
import time
import sys
import os
# import dask.dataframe as dd
import subprocess
import argparse
import re
from joblib import Parallel, delayed
import polars as pl


def get_args():
    """Function to get all args
    Arguments
    ---------

    
    meth_read = os.path.abspath(sys.argv[1])
    output_directory = os.path.abspath(sys.argv[2])
    genome_path = os.path.abspath(sys.argv[3])
    """
    parser = argparse.ArgumentParser(description="Generate bam file from ONT reads data",
    formatter_class=argparse.MetavarTypeHelpFormatter)
    parser.add_argument("--file","-f", help="Path to file with reads informations for one chromosome.", type=str, required=True, metavar="chr_file")
    parser.add_argument("--genome","-g", help="Path to the genome fasta file (Default : Col-CEN)", type=str, metavar="genome", required = False, default = "/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta")
    parser.add_argument("--basename_output","-base", help="Output file without extension", type=str, required = True)
    parser.add_argument("--output","-o", help="Path to output.", type=str, default = ".", required = False)
    parser.add_argument("--fast", help="Loading data into memory to work faster (Relatively high memory usage, need at least 2.5*FILE_SIZE Go of RAM).", default=False, action="store_true")
    parser.add_argument("--all", help="Compute all chromosomes in parallel.", default=False, action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        sys.exit("File doesn't exist")

    if not os.path.exists(args.genome):
        sys.exit("Genome file doesn't exist")

    if args.all:
        if "chr1" not in args.basename_output:
            sys.exit("Make sure that there is 'chr1' in the output basename")

    return args



def read_ref(reference_path, chrom):
    genome = {}
    with open(reference_path) as filin:
        for line in filin:
            if line.startswith(">"):
                chr = line[1:].strip()
                if chr == chrom:
                    genome[chr] = []
            else:
                if chr == chrom:
                    genome[chr] += list(line.strip())
    return genome

def get_chrs_path(args, sam_output, bam_basename):
    all_files = []
    path = os.path.abspath(args.file)
    print("Processing all chromosomes in parallel.")
    for i in range(1,6):
        bedfile = re.sub("chr\d",f"chr{i}",path)
        sam = re.sub("chr\d",f"chr{i}",sam_output)
        bam = re.sub("chr\d",f"chr{i}",bam_basename)

        if not os.path.exists(bedfile):
            print(f"{path} does not exists")
        else:
            all_files.append([bedfile, sam, bam])
    return all_files



def get_time(start):
    return round(time.time() - start,2)

def get_dict(read_file, chrom_name):
    header = ["chr","start","start_in_read","read","strand","state"]
    start_time = time.time()
    print(f"{chrom_name} :: Loading dataframe...")
    start_pl = time.time()
    meth_dtf = pl.read_csv(read_file, new_columns  = header, sep = "\t", has_header = False, low_memory = True , n_threads = 10).drop("start_in_read").to_pandas()
    start_opti = time.time()
    print(f"{chrom_name} :: Dataframe loaded in {get_time(start_pl)}s") 
    print(f"{chrom_name} :: Start memory optimisation") 
    meth_dtf["chr"] = meth_dtf["chr"].apply(lambda chrom: chrom[-1]).astype(np.int8)
    meth_dtf["strand"] = meth_dtf["strand"].astype("category")
    meth_dtf["state"] = meth_dtf["state"].astype(np.int8)
    print(f"{chrom_name} :: Memory optimisation done in {get_time(start_opti)}s") 
    print(f"{chrom_name} :: Converting dataframe to dictionary")
    start_time = time.time()
    df_dict = meth_dtf.set_index('read').groupby('read').apply(lambda x : x.to_numpy()).to_dict()
    print(f"{chrom_name} :: Dataframe converted into dictionary in {get_time(start_time)}s")
    return df_dict

def change_base(sequence, positions, base):
    for i in positions:
        if sequence[i] in ["T","A"]:
            sequence[i] = "N"
        else:
            print("ok")
            sequence[i] = base
    return sequence

def write_sam(read, flag, chrom, start, seq, file_out):
    before_seq = '*\t0\t0'
    qual = 5
    len_seq = len(seq)
    cigar = f"{len_seq}M"
    seq_quality = "#"*len_seq
    new_line = f"{read}\t{flag}\t{chrom}\t{start + 1}\t{qual}\t{cigar}\t{before_seq}\t{seq}\t{seq_quality}\n"
    file_out.write(new_line)

    
def create_sam_fast(df_dict, SAM_output, BAM_output, genome, args, chrom_name):
    bases_strand = {"+":"T", "-":"A"}
    flag_strand = {"+":  1, "-" : 17}
    with open(SAM_output, "w") as sam_filout:
        for read in df_dict:
            read_lines = df_dict[read]
            state = read_lines[:,3]
            base_to_change = np.where(state == 0)[0]
            pos = read_lines[:,1]
            first_position = pos.min()
            last_position = pos.max() + 1
            position_in_seq = (pos - first_position)[base_to_change]
            strand = read_lines[0,2]
            chrom = f"Chr{read_lines[0,0]}"
            base = bases_strand[strand]
            sequence = genome[chrom][first_position:last_position]
            sequence = "".join(change_base(sequence, position_in_seq, base))
            write_sam(read, flag_strand[strand], chrom, first_position, sequence, sam_filout)
    create_bam(BAM_output, SAM_output, args, chrom_name)



def create_sam_slow(bedfile, SAM_output, BAM_output, genome, args, chrom_name):
    
    # SORT IF NEEDED
    prev_read = ""
    act_read = ""
    bases_strand = {"+":"T", "-":"A"}
    flag_strand = {"+":  1, "-" : 17}
    pos_to_change = np.array([], dtype = np.int8)

    with open(bedfile) as bed_filin, open(SAM_output, "w") as sam_filout:
        lines = bed_filin.readline()
        items_first = lines.split()
        prev_read = items_first[3]
        prev_state = items_first[-1]
        first_position = int(items_first[1]) # assume to have bedfile sorted k1 k4 k2

        if prev_state == "0":
            pos_to_change = np.append(pos_to_change, 0)
            
        for line in bed_filin:
            items = line.split()
            act_read = items[3]
            if prev_read != act_read:
                last_position = prev_start + 1
                base = bases_strand[prev_strand]
                seq = genome[prev_chrom][first_position:last_position]
                seq = "".join(change_base(seq, pos_to_change, base))
                write_sam(prev_read, flag_strand[prev_strand], prev_chrom, first_position, seq, sam_filout)
                
                pos_to_change = np.array([], dtype = np.int8)
                prev_read = act_read
                first_position = int(items[1])
                
            # prev line info
            prev_chrom = items[0]
            prev_strand = items[-2]
            prev_state = items[-1]
            prev_start = int(items[1])
            if prev_state == "0":
                pos_to_change = np.append(pos_to_change, prev_start - first_position)
    

        # Process last line
        last_position = prev_start + 1
        base = bases_strand[prev_strand]
        seq = genome[prev_chrom][first_position:last_position]
        seq = "".join(change_base(seq, pos_to_change, base))
        write_sam(prev_read, flag_strand[prev_strand], prev_chrom, first_position, seq, sam_filout)
    create_bam(BAM_output, SAM_output, args, chrom_name)



def create_bam(bam_output, output_sam, args, chrom_name):

    genome_basename =  os.path.basename(args.genome)
    genome_name = os.path.splitext(genome_basename)[0]
    genome_dir = os.path.dirname(output_sam)
    chrom_sizes = f"{genome_dir}/{genome_name}_chrsize.txt"
    if not os.path.exists(chrom_sizes):
        cmd_chromsize_samtools = f"samtools faidx {args.genome}"
        print(f"{chrom_name} :: {cmd_chromsize_samtools}")
        subprocess.run([cmd_chromsize_samtools], shell=True)
        cmd_chromsize_cut = f"cut -f1,2 {args.genome}.fai > {chrom_sizes}"
        print(f"{chrom_name} :: {cmd_chromsize_cut}")
        subprocess.run([cmd_chromsize_cut], shell=True)

    header_bam = f"{os.path.dirname(bam_output)}/header.bam"
    with open(chrom_sizes) as filin, open(header_bam, "w") as filout:
        for line in filin:
            items = line.split()
            chrom = items[0]
            size = items[1]
            line_header = f"@SQ\tSN:{chrom}\tLN:{size}\n" # header bam
            filout.write(line_header)
    paste_header = f"cat {header_bam} {output_sam} > temp; mv temp {output_sam}"
    print(f"{chrom_name} :: {paste_header}")
    subprocess.run([paste_header], shell=True)

    samtools_bam = f"samtools view -b {output_sam} > {bam_output}"
    print(f"{chrom_name} :: {samtools_bam}")
    subprocess.run([samtools_bam], shell=True)

    samtools_sort = f"samtools sort {bam_output} > {bam_output[:-4]}.sort.bam"
    print(f"{chrom_name} :: {samtools_sort}")
    subprocess.run([samtools_sort], shell=True)

    samtools_index = f"samtools index {bam_output[:-4]}.sort.bam"
    print(f"{chrom_name} :: {samtools_index}")
    subprocess.run([samtools_index], shell=True)


def main_fast(bedfile_path, SAM_output, BAM_output, genome_path, args):
    begin_chr = time.time()
    try:
        chrom_name = re.search("chr\d", bedfile_path).group(0).capitalize()
    except:
        chrom_name = ""
    genome = read_ref(genome_path, chrom_name)
    df_dict = get_dict(bedfile_path, chrom_name)
    print(f"{chrom_name} :: Creating SAM file...")
    create_sam_fast(df_dict, SAM_output, BAM_output, genome, args, chrom_name)
    print(f"{chrom_name} :: Done for {chrom_name} in {get_time(begin_chr)}")


def main_slow(bedfile_path, SAM_output, BAM_output, genome_path, args):
    begin_chr = time.time()
    try:
        chrom_name = re.search("chr\d", bedfile_path).group(0).capitalize()
    except:
        chrom_name = ""
    genome = read_ref(genome_path, chrom_name)
    create_sam_slow(bedfile_path, SAM_output, BAM_output, genome, args, chrom_name)
    print(f"{chrom_name} :: Done for {chrom_name} in {get_time(begin_chr)}")


if __name__ == "__main__":
    # python /mnt/data2/rradjas/scripts/ont_to_bam.py --file meth_chr1.bed --genome /mnt/data2/rradjas/genome/Col-CEN/Col-CEN_v1.2.fasta --bam_output Cmt3_chr1.bam --output .
    start_time = time.time()
    
    args = get_args()
    bedfile = args.file
    basename = args.basename_output
    bam_folder = f"{args.output}/biseq_igv/bam"
    BAM_output = f"{args.output}/biseq_igv/bam/{basename}.bam"
    SAM_output = f"{args.output}/biseq_igv/{basename}.sam"
    os.makedirs(bam_folder, exist_ok=True)

    if args.fast:
        print("Running fast mode")
        run_function = main_fast
    else:        
        print("Running slow mode")
        run_function = main_slow


    if args.all:
        list_files = get_chrs_path(args, SAM_output, BAM_output)
        N_cpu = len(list_files)
        print(f"Using {N_cpu} CPUs")
        Parallel(n_jobs = N_cpu, verbose = 0, prefer="processes")(delayed(run_function)
        (bedfile_path, SAM_path, BAM_path, args.genome, args) for bedfile_path, SAM_path, BAM_path in list_files)
    else:
        run_function(bedfile, SAM_output, BAM_output, args.genome, args)

    print(time.time() - start_time)