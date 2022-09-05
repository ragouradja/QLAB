import re
import pandas as pd
import time
import numpy as np
import polars as pl
import sys


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
    meth_file = sys.argv[1] # meth_chr1.bed
    output = sys.argv[2]
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