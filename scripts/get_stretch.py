"""
Only one chr in input 
#Tombo:
reference genome : split by line ?


Install polars and pyarrow

Dont forget to sort the read file ! Otherwise the - strand will be reversed

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from scipy import stats
import polars as pl
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

    parser.add_argument("--read_file","-r", help="Path to the file containing read informations.", type=str, required=True, metavar="")
    parser.add_argument("--meth_file","-m", help="Path to the file containing methylation frequencies.", type=str, required=True, metavar="")
    parser.add_argument("--output_file","-o", help="Output file name", type=str, required = True, metavar="")
    args = parser.parse_args()
    return args




def get_dict_read(read_file):
    print(read_file)
    header = ["chr","start","end","read","strand","state"]
    start_time = time.time()
    dtf = pl.read_csv(read_file, new_columns  = header, sep = "\t", has_header = False, n_threads = 30).to_pandas()
    print("Read dtf loaded in : ",time.time() - start_time)

    start_dict = time.time()
    dict_dtf = dtf.groupby('read').apply(lambda g: list(map(list,g.values))).to_dict()
    print("Read dtf dict : ",time.time() - start_dict)

    # start_dict = time.time()
    # meth = dtf[["start","state"]]
    # dict_meth = meth.set_index("start").groupby("start").apply(lambda g: g.values).to_dict()
    # print("Meth dtf dict : ",time.time() - start_dict)
    
    return dict_dtf

def get_dict_pvalue(read_file):
    header = ["chr","start","end","strand", "frac","cov"]
    start_time = time.time()
    dtf = pl.read_csv(read_file, new_columns  = header, sep = "\t", has_header = False, n_threads = 30).to_pandas()
    print("Meth dtf loaded in : ",time.time() - start_time)

    start_pl = time.time()
    # print("Meth dtf to pandas : ",time.time() - start_pl)
    start_dict = time.time()
    meth = dtf[["start","frac","cov"]]
    dict_meth = meth[["start","frac","cov"]].set_index("start").to_dict(orient = "index")
    print("Meth dtf dict : ",time.time() - start_dict)

    return dict_meth


if __name__ == "__main__":
    args = get_args()
    read_file = args.read_file
    meth_bedgraph_file = args.meth_file # /mnt/data5/rradjas/ONT/Col-0/bedgraph --> chr start end frac cov
    output_file = args.output_file
    print(f"Writing in {output_file}")
    start_time = time.time()


    dict_stretch =  get_dict_read(read_file)
    dict_meth = get_dict_pvalue(meth_bedgraph_file)

    with open(output_file,"w") as f:
        for read in dict_stretch:
            read_stretch = dict_stretch[read]
            prev_zero = True
            for line in read_stretch:
                state = line[5]
                if state == 0 or line == read_stretch[-1]: # end of stretch or no stretch or end of read
                    if not prev_zero and length >= 5: # end of stretch <=> after "1"
                        cov /= length
                        pvalue = round(-np.log10(stats.binom_test(1, n=round(cov), p=frac, alternative='greater')),2)
                        content = f"{line[0]}\t{stretch_start}\t{stretch_end}\t{read}\t{line[4]}\t{length}\t{pvalue}\n"
                        f.write(content)
                    prev_zero = True
                else: # state == 1
                    if prev_zero:
                        stretch_start = line[1]
                        length = 1
                        cov = dict_meth[stretch_start]["cov"]
                        frac = dict_meth[stretch_start]["frac"]
                        prev_zero = False
                    else:
                        length += 1
                        cov += dict_meth[line[1]]["cov"]
                        frac *= dict_meth[line[1]]["frac"]
                    stretch_end = line[1]
                    
    print("Total time : ",time.time() - start_time)
