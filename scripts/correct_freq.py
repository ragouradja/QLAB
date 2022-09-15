"""Getting methylation frequencies"""

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
    parser = argparse.ArgumentParser(description="Getting frequencies of methylation",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--file", "-f", help="Path to the bedfile with reads information", required = True, type = str)
    parser.add_argument("--output_file", "-o", help="Path to output file", required = True, type = str, default=".")
    parser.add_argument("--one_chr", "-all", help="If the file contains only one chromosome (faster)", required = False,  default= False, action="store_true")
    args = parser.parse_args()  
    
    if not os.path.exists(args.file):
        sys.exit("File doesn't exist")
    return args

def get_dict(meth_file):
    header = ["chr","start","end","read","strand","state"]
    start_time = time.time()

    dtf = pl.read_csv(meth_file, new_columns  = header, sep = "\t", has_header = False, n_threads = 40).drop(["end","read"]).to_pandas()
    print("load",time.time() - start_time)

    start_time = time.time()
    if args.one_chr:
        dtf["id"] = dtf['chr'] + "_" + dtf['start'].astype(str)
        dict_meth = dtf.set_index("id").groupby("id").apply(lambda g: g.values).to_dict()
    else:
        dict_meth = dtf.set_index("start").groupby("start").apply(lambda g: g.values).to_dict()
    
    print("dict",time.time() - start_time)
    return dict_meth

def freq_all_chr(dict_meth, output):
    with open(output, "w") as fout:
        for chr_pos in dict_meth.keys():
            lines = dict_meth[chr_pos]
            chrom = lines[0][0]
            pos = lines[0][1]
            strand = lines[0][-2]
            cov = lines.shape[0]
            meth = round(lines[:,-1].sum() / cov,4)
            content = f"{chrom}\t{pos}\t{pos + 1}\t{strand}\t{meth}\t{cov}\n"    
            fout.write(content)

def freq(dict_meth, output):
    with open(output, "w") as fout:
        for pos in dict_meth.keys():
            lines = dict_meth[pos]
            chrom = lines[0][0]
            strand = lines[0][-2]
            cov = lines.shape[0]
            meth = round(lines[:,-1].sum() / cov,4)
            content = f"{chrom}\t{pos}\t{pos + 1}\t{strand}\t{meth}\t{cov}\n"    
            fout.write(content)

if __name__ == "__main__":

    start_time = time.time()    
    args = get_args()
    meth_file = args.file
    output = args.output_file
    dict_meth = get_dict(meth_file)   

    if args.one_chr:
        freq(dict_meth, output)
    else:
        freq_all_chr(dict_meth, output)

    print("Done in ", time.time() - start_time)