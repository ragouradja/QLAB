"""Separate TSV file in chromosomes and in context (CG, CHG, CHH, A)"""

import sys
import argparse
import os
import re
import time
import gzip


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

    parser.add_argument("--file", "-f", help="Path to BED or TSV file to separate", required = True, type = str, metavar = "")
    parser.add_argument("--genome", "-g", help="Path to genome file", required = False, type = str, metavar = "", default="/mnt/data2/rradjas/genome/Col-CEN/fasta/Col-CEN_all.fasta")
    parser.add_argument("--output_dir", "-o", help="Output directory to create folders for each chromosom", required = True, default=".", type = str, metavar = "")
    parser.add_argument("--chrom_output", help="Name of the output files for chr1 file", required = False, type = str, metavar = "", default= "chr1_methylation.bed")
    parser.add_argument("--context_output", help="Name of the output files for CG file", required = False, type = str, metavar = "", default= "CG_context.bed")
    parser.add_argument('--awk', nargs='+', help='Index of columns to keep (1-based index)', required=False, default=None, type=int)
    parser.add_argument("--strand_column", "-s", help="Index of column containing strand information (1-based index)", required = False, default=5, type = int, metavar = "")
    parser.add_argument("--tsv", "-t", help="Processing tsv file from DSP pipeline", required = False, action="store_true")

    args = parser.parse_args()

    if not os.path.exists(args.file):
        sys.exit(f"{args.file} doesn't exist")
    if not os.path.exists(args.genome):
        sys.exit(f"{args.genome} doesn't exist")
    if args.context_output != "CG_context.bed":
        if "CG" not in args.context_output:
            sys.exit("Please make sure that there is 'CG' in your context output filename")
    if args.chrom_output != "chr1_methylation.bed":
        if "chr1" not in args.chrom_output:
            sys.exit("Please make sure that there is 'chr1' in your chromosom output filename")

    for i in range(1,6):
        chrom = f"Chr{i}"
        os.makedirs(f"{args.output_dir}/{chrom.lower()}/context", exist_ok=True)
    return args


def read_ref(reference_path):
    genome = {}
    print("Reading genome...")
    with open(reference_path) as filin:
        for line in filin:
            if line.startswith(">"):
                chrom = line[1:].strip()
                genome[chrom] = ""
            else:
                genome[chrom] += line.strip()
    print("Done")
    return genome

def get_context(chrom, pos, strand, genome):
    H = ["A","T","C"]
    Hneg = ["A","T","G"]
    if strand == "+":
        kmer = get_kmer_pos(chrom, pos, genome)
        if kmer[0] == "A":
            return "A"
        if kmer[1] in H:
            if kmer[2] in H:
                return "CHH"
            return "CHG"
        return "CG"
    elif strand == "-":
        kmer = get_kmer_neg(chrom, pos, genome)
        if kmer[2] == "T":
            return "A"
        if kmer[1] in Hneg:
            if kmer[0] in Hneg:
                return "CHH"
            return "CHG"
        return "CG"
    else:
        sys.exit("Please specify the right column position with the strand information")

    
def get_kmer_pos(chrom, pos, genome):
    first = pos 
    last = first + 3
    return genome[chrom][first:last]

def get_kmer_neg(chrom, pos, genome):
    first = pos - 2
    last = pos + 1
    return genome[chrom][first:last]

def read_bed(genome, args):
    all_files_5mc = {}
    all_files_chrom = {}
    for i in range(1,6):
        chrom = f"Chr{i}"
        chrom_file = re.sub("chr1", chrom.lower(), args.chrom_output)
        all_files_5mc[chrom] = {}
        all_files_chrom[chrom] = open(f"{args.output_dir}/{chrom.lower()}/{chrom_file}", "w")

    file = args.file
    if file.endswith(".gz"):
        filin = gzip.open(file, 'rt')
    else:
        filin = open(file, "r")

    for line in filin:
        items = line.split()
        chrom = items[0]
        start = int(items[1])
        strand = items[args.strand_column - 1] # from 1-based --> 0-based index

        if chrom in ["ChrC", "ChrM"]:
            continue

        if args.awk:
            line = ""
            for position in args.awk:
                line += f"{items[position-1]}\t"
            line = line[:-2] + "\n"
        elif args.tsv:
            line = f"{chrom}\t{start}\t{start+1}\t{items[4]}\t{strand}\t{items[8]}\n"
        context = get_context(chrom, start, strand, genome)
        if context not in all_files_5mc[chrom]:
            out_file_context = re.sub("CG", context, args.context_output)
            all_files_5mc[chrom][context] = open(f"{args.output_dir}/{chrom.lower()}/context/{out_file_context}", "w")

        all_files_5mc[chrom][context].write(line)
        all_files_chrom[chrom].write(line)

    filin.close()
    for chrom in all_files_chrom:
        all_files_chrom[chrom].close()
        for context in all_files_5mc[chrom]:
            all_files_5mc[chrom][context].close()


if __name__ == "__main__":
    start_time = time.time()
    args = get_args()
    genome_ref = read_ref(args.genome)
    read_bed(genome_ref, args)
    end_time = time.time() - start_time
    print(f"Done in {end_time}s")    