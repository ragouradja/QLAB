""" Converting BED file to BAM file with Biseq mode visualization """



import pandas as pd
import numpy as np
import time
import sys
import os
import subprocess
import argparse
import re
from joblib import Parallel, delayed
import polars as pl


def get_args():
    """Function to get all args"""

    parser = argparse.ArgumentParser(description="Generate bam file from ONT reads data",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--file","-f", help="Path to file with reads informations for one chromosome.", type=str, required=True, metavar="chr_file")
    parser.add_argument("--genome","-g", help="Path to the genome fasta file (Default : Col-CEN)", type=str, metavar="genome", required = False, default = "/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta")
    parser.add_argument("--basename_output","-base", help="Output file without extension", type=str, required = True)
    parser.add_argument("--output","-o", help="Path to output directory.", type=str, default = ".", required = False)
    parser.add_argument("--fast", help="Loading data into memory to work faster (Relatively high memory usage, need at least 2.5*FILE_SIZE Go of RAM).", default=False, action="store_true")
    parser.add_argument("--all", help="Compute all chromosomes in parallel.", default=False, action="store_true")
    parser.add_argument("--all_in_one", help="If your file has all chromosomes.", default=False, action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        sys.exit("File doesn't exist")

    if not os.path.exists(args.genome):
        sys.exit("Genome file doesn't exist")
    return args



def read_ref(reference_path, chrom, args):
    """Read the genome
    Arguments
    ---------
    reference_path : string
        Path to the reference genome
    chrom : List
        List of chromosom to keep while reading the genome 
        (can save some memory space if some chrom are not used)
    args : List
        List of arguments from get_args()
        

    Return
    ------
        genome : Dict
            Genome in format : {"Chr1" : ["A","T","C"...],
                                "Chr2" : [...]}
    """

    genome = {}
    print("Reading genome...")
    with open(reference_path) as filin:
        for line in filin:
            if line.startswith(">"):
                chr = line[1:].strip()
                if args.all_in_one or chr == chrom:
                    genome[chr] = []
            else:
                if args.all_in_one or chr == chrom:
                    genome[chr] += list(line.strip())
    print("Done")
    return genome

def get_chrs_path(args):
    """Find path to all chromosomes bedfiles from chr1 bedfile path
    using regular expression to replace 'chr1' with 'chr2','chr3' ...

    This function is used only if you run with --all option.

    Arguments
    ---------
    args : List
        List of arguments from get_args()
        

    Return
    ------
        all_files : List
            List containing existing paths to bedfiles for each chrom
    """
    all_files = []
    path = os.path.abspath(args.file)
    print("Processing all chromosomes in parallel.")
    for i in range(1,6):
        bedfile = re.sub("chr\d",f"chr{i}",path)

        if not os.path.exists(bedfile):
            print(f"{path} does not exists")
        else:
            all_files.append(bedfile)
    return all_files



def get_time(start):
    return round(time.time() - start,2)

def get_dict(read_file, chrom_name):
    """Read methylation data (for fast_mode only)
    
    Reading data using Polars which load data using multiprocessing.
    Optimization of memory by using cheaper data type.
    Groupby data by read then converting the dataframe into dictionary which use less memory.

    Arguments
    ---------
    read_file : string
        Path to the bedfile containing methylation information
        per read and per position.
    chrom_name : string
        Chromosome name (only for printing informations)

    Return
    ------
        df_dict : Dict
            Methylation data in dict format : {Read1 : np.array(["chr","start","strand","state"]), 
                                               Read2 : ...}
            Remark : Data for each read are in numpy array which allow vectorial operations
                     within one read (faster than classical operations).
    """

    header = ["chr","start","start_in_read","read","strand","state"]
    start_time = time.time()
    print(f"{chrom_name} :: Loading dataframe...")
    start_pl = time.time()
    meth_dtf = pl.read_csv(read_file, new_columns  = header, sep = "\t", has_header = False,
     low_memory = True , n_threads = 10).drop("start_in_read").to_pandas()
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
    """Converting sequence to Biseq mode
    
    Arguments
    ---------
    sequence : List
        Sequence to modify

    positions : List
        Position of unmethylated cytosines to modify
    base : string
        New base to replace cytosine
        ('T' if + strand, 'A' if - strand)

    Return
    ------
    sequence : List
        Sequence modified 
    """

    for i in positions:
        sequence[i] = base
    return sequence

def write_sam(read, flag, chrom, start, seq, file_out):
    """ Function to write one line in SAM file"""

    before_seq = '*\t0\t0'
    qual = 5
    len_seq = len(seq)
    cigar = f"{len_seq}M"
    seq_quality = "#"*len_seq
    new_line = f"{read}\t{flag}\t{chrom}\t{start + 1}\t{qual}\t{cigar}\t{before_seq}\t{seq}\t{seq_quality}\n"
    file_out.write(new_line)

    
def create_sam_fast(df_dict, SAM_output, BAM_output, genome, args, chrom_name):
    """Creating SAM file in Biseq mode (fast_mode and relatively high memory usage)
    
    Iterate through all keys of df_dict, which are all reads, and modify unmethylated cytosines.
    The sequences modifed are not exactly the raw sequence of the read 
        but directly extracted from the genome (to avoid dealing with insertion or deletion)

    Arguments
    ---------
    df_dict : dict
        Methylation data in dict format from get_dict()
    SAM_output : string
        Filename of SAM file
    BAM_output : string
        Filename of BAM file
    genome : dict
        Genome in dict format from read_ref()
    args : list
        List of arguments from get_args()
        
    chrom_name : string 
        Chromosome name (only for printing informations)


    """


    bases_strand = {"+":"T", "-":"A"}
    flag_strand = {"+":  1, "-" : 17}
    with open(SAM_output, "w") as sam_filout:
        for read in df_dict:
            read_lines = df_dict[read]
            # Get methylation information for this particular read
            state = read_lines[:,3]
            # Selecting unmethylated cytosines
            base_to_change = np.where(state == 0)[0]
            pos = read_lines[:,1]
            first_position = pos.min()
            last_position = pos.max() + 1
            # Converting relative position to absolute position (starting at 0)
            position_in_seq = (pos - first_position)[base_to_change]
            strand = read_lines[0,2]
            chrom = f"Chr{read_lines[0,0]}"
            # Base that will replace cytosines, depending on strand
            base = bases_strand[strand]
            # Extracting sequence from genome dict
            # The first (also for last) nucleotide in the sequence 
            # will be the first (last) cytosine detected for a given read
            sequence = genome[chrom][first_position:last_position]
            # Modifying sequence
            sequence = "".join(change_base(sequence, position_in_seq, base))
            # Write one line in SAM file for this read with modified sequence
            write_sam(read, flag_strand[strand], chrom, first_position, sequence, sam_filout)



def create_sam_slow(bedfile, SAM_output, BAM_output, genome, args, chrom_name):
    """Creating SAM file in Biseq mode (slow_mode and no memory usage)
    Need to have bedfiles sorted by k4 k2 : By read and then by the position

    Iterate through all lines of one chromosome file and modify the sequence.
    Since we are not loading all the data into memory, we can parallelize
        the reading and processing of all 5 bedfiles.
    The sequences modifed are not exactly the raw sequence of the read 
        but directly extracted from the genome (to avoid dealing with insertion or deletion)

    Arguments
    ---------
    bedfile : string
        Path to one bedfile to read and process
    SAM_output : string
        Filename of SAM file
    BAM_output : string
        Filename of BAM file
    genome : dict
        Genome in dict format from read_ref()
    args : list
        List of arguments from get_args()
    chrom_name : string 
        Chromosome name (only for printing informations)


    """
    prev_read = ""
    act_read = ""
    bases_strand = {"+":"T", "-":"A"}
    flag_strand = {"+":  1, "-" : 17}
    pos_to_change = np.array([], dtype = np.int8)

    # Same logic as for create_sam_fast() but here we iterate directly the bedfile
    # To detect the end of one read, we need to keep in memory the last line ('prev' variables)
    with open(bedfile) as bed_filin, open(SAM_output, "a") as sam_filout:
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
            # If we detect a new read, the sequence of the previous read is modified
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
    



def create_bam(bam_output, output_sam, args):
    """Write BAM file from SAM file (from fast or slow mode)"""

    genome_basename =  os.path.basename(args.genome)
    genome_name = os.path.splitext(genome_basename)[0]
    genome_dir = os.path.dirname(output_sam)
    chrom_sizes = f"{genome_dir}/{genome_name}_chrsize.txt"
    if not os.path.exists(chrom_sizes):
        cmd_chromsize_samtools = f"samtools faidx {args.genome}"
        print(f"{cmd_chromsize_samtools}")
        subprocess.run([cmd_chromsize_samtools], shell=True)
        cmd_chromsize_cut = f"cut -f1,2 {args.genome}.fai > {chrom_sizes}"
        print(f"{cmd_chromsize_cut}")
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
    print(f"{paste_header}")
    subprocess.run([paste_header], shell=True)

    samtools_bam = f"samtools view -b {output_sam} > {bam_output}"
    print(f"{samtools_bam}")
    subprocess.run([samtools_bam], shell=True)

    samtools_sort = f"samtools sort {bam_output} > {bam_output[:-4]}.sort.bam"
    print(f"{samtools_sort}")
    subprocess.run([samtools_sort], shell=True)

    samtools_index = f"samtools index {bam_output[:-4]}.sort.bam"
    print(f"{samtools_index}")
    subprocess.run([samtools_index], shell=True)


def main_fast(bedfile_path, SAM_output, BAM_output, genome_path, args):
    begin_chr = time.time()
    try:
        chrom_name = re.search("chr\d", bedfile_path).group(0).capitalize()
    except:
        chrom_name = ""
    genome = read_ref(genome_path, chrom_name, args)
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
    genome = read_ref(genome_path, chrom_name, args)
    create_sam_slow(bedfile_path, SAM_output, BAM_output, genome, args, chrom_name)
    print(f"{chrom_name} :: Done for {chrom_name} in {get_time(begin_chr)}")


if __name__ == "__main__":
    start_time = time.time()
    
    args = get_args()
    bedfile = args.file
    basename = args.basename_output
    bam_folder = f"{args.output}/biseq_igv/bam"
    BAM_output = f"{args.output}/biseq_igv/bam/{basename}.bam"
    SAM_output = f"{args.output}/biseq_igv/{basename}.sam"
    if os.path.exists(SAM_output):
        os.remove(SAM_output)
    os.makedirs(bam_folder, exist_ok=True)

    if args.fast:
        print("Running fast mode")
        run_function = main_fast
    else:        
        print("Running slow mode")
        run_function = main_slow


    if args.all:
        list_files = get_chrs_path(args)
        N_cpu = len(list_files)
        print(f"Using {N_cpu} CPUs")
        Parallel(n_jobs = N_cpu, verbose = 0, prefer="processes")(delayed(run_function)
        (bedfile_path, SAM_output, BAM_output, args.genome, args) for bedfile_path in list_files)
    else:
        run_function(bedfile, SAM_output, BAM_output, args.genome, args)

    create_bam(BAM_output, SAM_output, args)

    print(time.time() - start_time)