"""Separate methylation data into 4 context CG CHG, CHH and A from genome reference and bedfile
"""
import sys
import argparse
import os
import time
import re

def get_args():
    """Function to get all args
    Arguments
    ---------

    
    meth_read = os.path.abspath(sys.argv[1])
    output_directory = os.path.abspath(sys.argv[2])
    genome_path = os.path.abspath(sys.argv[3])
    """
    parser = argparse.ArgumentParser(description="Separate file according to chromosom or/and methylation context",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--file", "-f", help="Path to methylation file", required = True, type = str, metavar = "")
    parser.add_argument("--genome", "-g", help="Path to genome file", required = False, type = str, metavar = "", default="/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta")
    parser.add_argument("--output_dir", "-o", help="Output directory to create file for each context", required = True, type = str, metavar = "")
    parser.add_argument("--strand_column", "-s", help="Index (1-based) of column containing strand information", required = False, default=5, type = int, metavar = "")
    parser.add_argument("--output_file", "-of", help="Output filename example for CG context (e.g CG_context.bed)", required = False, type = str, metavar = "", default="CG_context.bed")

    args = parser.parse_args()

    if not os.path.exists(args.file):
        sys.exit(f"{args.file} doesn't exist")
    if not os.path.exists(args.genome):
        sys.exit(f"{args.genome} doesn't exist")
    if not os.path.exists(args.output_dir):
        print(f"Creating {args.output_dir} folder")
        os.makedirs(args.output_dir)
    return args




def read_ref(reference_path):
    genome = {}
    with open(reference_path) as filin:
        for line in filin:
            if line.startswith(">"):
                chr = line[1:].strip()
                genome[chr] = ""
            else:
                genome[chr] += line.strip()
    return genome

def get_context(chr, pos, strand, genome):
    H = ["A","T","C"]
    Hneg = ["A","T","G"]
    if strand == "+":
        kmer = get_kmer_pos(chr, pos, genome)
        if kmer[0] == "A":
            return "A"
        if kmer[1] in H:
            if kmer[2] in H:
                return "CHH"
            return "CHG"
        return "CG"
    elif strand == "-":
        kmer = get_kmer_neg(chr, pos, genome)
        if kmer[2] == "T":
            return "A"
        if kmer[1] in Hneg:
            if kmer[0] in Hneg:
                return "CHH"
            return "CHG"
        return "CG"
    else:
        sys.exit("Please specify the right column position with the strand information")


def get_kmer_pos(chr, pos, genome):
    first = pos 
    last = first + 3
    return genome[chr][first:last]

def get_kmer_neg(chr, pos, genome):
    first = pos - 2
    last = pos + 1
    return genome[chr][first:last]

def read_bed(args):
    start_time = time.time()
    path_bed = args.file
    path_ref = args.genome
    output_dir = args.output_dir
    strand_col = args.strand_column
    genome = read_ref(path_ref)
    mc5 = {}
    with open(path_bed) as bedfile:
        for line in bedfile:
            items = line.split()
            chr = items[0]
            strand = items[strand_col - 1] # from 1-based --> 0-based index
            pos = int(items[1])
            context = get_context(chr, pos, strand, genome)
            if context not in mc5:
                out_file = re.sub("CG", context, args.output_file)
                mc5[context] = open(f"{output_dir}/{out_file}", "w")
            mc5[context].write(line)

    for context in mc5:
        mc5[context].close()
    print("Total time : ", time.time() - start_time)

if __name__ == "__main__":
    args = get_args()
    read_bed(args)
