"""Bam to bed from sam file"""

import sys
import os
import time
import re
from joblib import Parallel, delayed
import argparse

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
    parser.add_argument("--output","-o", help="Path to output.", type=str, default = ".", required = False)
    args = parser.parse_args()

    if not os.path.exists(args.file):
        sys.exit("File doesn't exist")

    if not os.path.exists(args.genome):
        sys.exit("Genome file doesn't exist")
    return args




def read_ref(reference_path, chrom):
    genome = {}
    with open(reference_path) as filin:
        for line in filin:
            if line.startswith(">"):
                chr = line[1:].strip()
                if chr == chrom:
                    genome[chr] = ""
            else:
                if chr == chrom:
                    genome[chr] += line.strip()
    print(genome.keys())
    return genome

def process_chrom(samfile, chrom, letter, genome_path, output_dir):
    Chrom = f"Chr{chrom[-1]}"
    chrom = f"chr{chrom[-1]}"
    genome = read_ref(genome_path, Chrom)
    contexts = ["CG","CHG","CHH"]
    meth_dico = {"Z":["CG",1], "X":["CHG",1],"H":["CHH",1], "z":["CG",0], "x":["CHG",0], "h":["CHH",0]}
    strand_dico = {"C":"+",  "G":"-"}
    files = {}
    
    reg = re.compile("chr\d")
    samfile_chrom = re.sub("chr\d",chrom,samfile)

    chrom_start = time.time()

    for context in contexts:
        files[context] = open(f"{output_dir}/{chrom}/context/{context}_biseq.bed","a")

    print(f"Reading sam file : {samfile_chrom}")   
    with open(samfile_chrom) as filin:
        for line in filin:
            items = line.split()
            if len(items[3]) > 4:
                continue
            seq = items[-1][5:]
            start_read = int(items[1])
            for i in range(len(seq)):
                char = seq[i]
                if char in meth_dico:
                    context, state = meth_dico[char]
                    start = start_read+i-1
                    end = start_read+i
                    read = items[2] + letter
                    strand = strand_dico[genome[Chrom][start]]

                    # if chrom in chromosoms: no need to check since I already filtered out ChrC and ChrM during chromosom speparation
                    files[context].write(f"{Chrom}\t{start}\t{end}\t{read}\t{strand}\t{state}\n")

    for context in contexts:
        files[context].close()

    print(f"Sam file {samfile_chrom} done in {time.time() - chrom_start}")


if __name__ == "__main__":
    start_time = time.time()
    args = get_args()
    samfile = args.file
    genome_path = args.genome
    
    output_dir = "reads_analysis"
    chromosoms = ["chr1","chr2","chr3","chr4","chr5"]
    if "reverse" in samfile:
        letter = "r"            
    else:
        letter = "f"         

    # Creating Folder chrom
    for chrom in chromosoms:
        if not os.path.exists(f"{output_dir}/{chrom}/context"):
            os.makedirs(f"{output_dir}/{chrom}/context", exist_ok=True)

    Parallel(n_jobs = 5, verbose = 0, prefer="processes")(delayed(process_chrom)
    (samfile, chrom, letter, genome_path, output_dir) for chrom in chromosoms)

    print("Total time: ", time.time() - start_time)


