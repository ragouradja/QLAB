import re
import pandas as pd
import time
import numpy as np
import polars as pl
import sys
import argparse
import os


def get_args():
    """Function to get all args
    Arguments
    ---------

    
    meth_read = os.path.abspath(sys.argv[1])
    output_directory = os.path.abspath(sys.argv[2])
    genome_path = os.path.abspath(sys.argv[3])
    """
    parser = argparse.ArgumentParser(description="Compute average methylation level per position",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--file", "-f", help="Path to bedfile with read information", required = True, type = str, metavar = "")
    parser.add_argument("--output_file", "-o", help="Output filename", required = True, type = str, metavar = "")

    args = parser.parse_args()

    if not os.path.exists(args.file):
        sys.exit(f"{args.file} doesn't exist")
    return args



def get_dict(meth_file):
    header = ["chr","start","end","read","strand","state"]
    start_time = time.time()

    dtf = pl.read_csv(meth_file, new_columns  = header, sep = "\t", has_header = False, n_threads = 40).to_pandas()
    print("load",time.time() - start_time)

    start_time = time.time()
    dict_meth = dtf.set_index("start").groupby("start").apply(lambda g: g.values).to_dict()
    print("dict",time.time() - start_time)
    return dict_meth

if __name__ == "__main__":

    start_time = time.time()    
    args = get_args()
    meth_file = args.file
    output = args.output_file
    print(meth_file, output)
    dict_meth = get_dict(meth_file)   
    with open(output, "w") as fout:
        for pos in dict_meth.keys():
            lines = dict_meth[pos]
            chrom = lines[0][0]
            strand = lines[0][3]
            cov = lines.shape[0]
            meth = round(lines[:,-1].sum() / cov,4)
            content = f"{chrom}\t{pos}\t{pos + 1}\t{strand}\t{meth}\t{cov}\n"    
            fout.write(content)

    print("Done in ", time.time() - start_time)