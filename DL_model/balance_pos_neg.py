import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import sklearn
import time
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
    parser = argparse.ArgumentParser(description="Getting frequencies of methylation",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pos_file", "-p", help="Path to the positive dataset", required = True, type = str)
    parser.add_argument("--neg_file", "-n", help="Path to the negative dataset", required = True, type = str)
    parser.add_argument("--context", "-c", help="Methylation context (CG, CHG or CHH)", required = True, type = str)
    args = parser.parse_args()  
    
    if not os.path.exists(args.pos_file):
        sys.exit(f"{args.pos_file} doesn't exist")

    if not os.path.exists(args.neg_file):
        sys.exit(f"{args.neg_file} doesn't exist")

    if args.context not in ["CG","CHG","CHH"]:
        sys.exit("Select context among CG, CHG or CHH")
    return args



def plot_kmer(kmer, context, save_path = None, plotting = True):
    if "pos" in save_path:
        state = "positive"
    else:
        state = "negative"

    kmer = kmer.groupby("kmer").count()["chrom"].sort_values()
    plt.figure(figsize = (12,8))
    plt.hist(kmer)
    plt.xlabel("Occurrence for a given kmer", fontsize = 12)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)

    plt.title(f"{context} {state} kmer", fontsize = 15)
    if save_path:
        plt.savefig(save_path)
    if plotting:
        plt.plot()


if __name__ == "__main__":
    start_time = time.time()
    args = get_args()
    pos_path = args.pos_file
    neg_path = args.neg_file
    context = args.context
    output_dir , neg_output = os.path.split(neg_path)
    output_dir = f"{output_dir}/balanced"
    output_dir_image = f"{output_dir}/images"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir_image, exist_ok=True)

    neg_output = f"{neg_output[:-4]}.balanced.tsv"
    header = ["chrom","start","strand","pos","read","t", "kmer","col8","col9","col10","col11","state"]
    pos_CG_raw = pl.read_csv(pos_path, new_columns  = header, sep = "\t", has_header = False).to_pandas()
    neg_CG_raw = pl.read_csv(neg_path, new_columns  = header, sep = "\t", has_header = False).to_pandas()
    print("Data loaded")
    pos_CG = pos_CG_raw[["chrom","start","kmer"]]
    neg_CG = neg_CG_raw[["chrom","start","kmer"]]

    plot_kmer(pos_CG, context, save_path = f"{output_dir_image}/{context}_raw_kmer_pos.pdf", plotting = False)
    plot_kmer(neg_CG, context, save_path = f"{output_dir_image}/{context}_raw_kmer_neg.pdf", plotting = False)

    Kpos = pos_CG.kmer.values
    Kneg = neg_CG.kmer.values
    print("Start intersect")

    Kcomm = np.intersect1d(Kneg, Kpos, assume_unique=False)
    print("Start diff")

    Kdiff = np.setdiff1d(Kneg, Kpos)
    KNUM_pos = pos_CG.groupby("kmer").count()["chrom"].sort_values()

    Sneg_= np.array([])
    for k in Kcomm:
        kcount = KNUM_pos[k]
        neg_CG_k = neg_CG[neg_CG.kmer == k]
        Sneg_k_index = np.random.choice(neg_CG_k.index, kcount)
        Sneg_ = np.append(Sneg_, Sneg_k_index)
        
    N = pos_CG.shape[0] - Sneg_.shape[0]
    if N > 0:
        Kdiff_neg = neg_CG[neg_CG.kmer.isin(Kdiff)]
        fill_kmer_index = np.random.choice(Kdiff_neg.index, N)
    Sneg_ = np.append(Sneg_, fill_kmer_index)

    final_neg = sklearn.utils.shuffle(neg_CG_raw.iloc[Sneg_])
    final_neg.to_csv(f"{output_dir}/{neg_output}", index = False, sep = "\t", header = False)
    print(f"New negative dataset saved in {output_dir}/{neg_output}")
    plot_kmer(final_neg, f"{context} balanced", save_path = f"{output_dir_image}/{context}_raw_kmer_neg_balanced.pdf", plotting = False)
    print("Done in ", time.time() - start_time)
