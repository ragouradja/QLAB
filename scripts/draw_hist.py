import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import os
import argparse
import sys
import numpy as np
import re
import pandas as pd


def get_args():
    """Function to get all args
    Arguments
    ---------

    
    meth_read = os.path.abspath(sys.argv[1])
    output_directory = os.path.abspath(sys.argv[2])
    genome_path = os.path.abspath(sys.argv[3])
    """
    parser = argparse.ArgumentParser(description="Plotting histogram",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--file", "-f", help="Path to file containing path for CG stretch file for each mutant (1 mutant per line : NAME PATH)", required = True, type = str)
    parser.add_argument("--outdir", "-o", help="Path to output folder", required = False, type = str, default="images")
    parser.add_argument("--hist", "-hi", help="Draw histogram", required = False, action = "store_true", default = False)
    parser.add_argument("--box", "-b", help="Draw boxplot", required = False, action = "store_true", default = False)
    parser.add_argument("--violin", "-vi", help="Draw violin plot", required = False, action = "store_true", default = False)
    parser.add_argument("--hist_context", "-hc", help="Plot histogram per context", required = False, default=False, action = "store_true")
    parser.add_argument("--title_hist", "-th", help="Title of the histogram", required = False, type = str, metavar = "", default="Normalized number of significant stretches in ONT mutants")
    parser.add_argument("--title_violin", "-tv", help="Title of the violin plot", required = False, type = str, metavar = "", default="Length of statistically significant CG stretches")

    parser.add_argument("--ylab_hist", "-yh", help="Label of y axis", required = False, type = str, metavar = "", default="Normalized amount")
    parser.add_argument("--ylab_violin", "-yv", help="Label of y axis", required = False, type = str, metavar = "", default="Length in bp (log10)")
    parser.add_argument("--ylab_box", "-yb", help="Label of y axis", required = False, type = str, metavar = "", default="-log10 pvalue")

    parser.add_argument("--output_hist", "-outh", help="Output filename of histogram with extension (pdf, png...)", required = False, type = str, metavar = "", default = "hist_stretch.pdf")
    parser.add_argument("--output_violin", "-outv", help="Output filename of violin plot with extension (pdf, png...)", required = False, type = str, metavar = "", default = "CG_stretch_violin.pdf")

    parser.add_argument("--xsize", "-xs", help="Size of xaxis label", required = False, type = int, metavar = "", default=15)
    parser.add_argument("--ysize", "-ys", help="Size of yaxis label", required = False, type = int, metavar = "", default=15)
    parser.add_argument("--yaxis_size", "-yts", help="Size of yaxis values", required = False, type = int, metavar = "", default=15)
    parser.add_argument("--legend_size", "-ls", help="Size of legend", required = False, type = int, metavar = "", default=15)
    

    args = parser.parse_args()  
    
    if args.hist_context:
        if "CG" not in args.output_hist:
            sys.exit("Make sure that 'CG' is in the output filename of hist plot")
    if "CG" not in args.output_violin:
        sys.exit("Make sure that 'CG' is in the output filename of violin plot")
    if not (args.hist or args.violin):
        sys.exit("Select --hist or/and --violin to create a plot")
    if args.outdir != ".":
        os.makedirs(args.outdir, exist_ok=True)
        args.output_hist = f"{args.outdir}/{args.output_hist}"
        args.output_violin = f"{args.outdir}/{args.output_violin}"

    return args

def get_data(args):
    dict_data = {}
    with open(args.file) as filin:
        for line in filin:
            if line.startswith("#"):
                continue
            items = line.split()
            if len(items) == 0:
                continue
            name = items[0]
            path = items[1]
            dict_data[name] = []
            if len(items) == 3:
                div = float(items[2])
            else:
                div = 1
            for context in ["CG","CHG","CHH"]:
                context_path = re.sub("CG",context, path)
                if not os.path.exists(context_path):
                    print(f"{context_path} doesn't exists. Skipping {context} context")
                    continue
                f = open(context_path)
                content = f.readlines()
                N_stretch = len(content) / div
                f.close()
                dict_data[name].append(N_stretch)
    return dict_data


def bar_plot(ax, data, args, colors=None, total_width=0.8, single_width=1, legend=True):
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    n_bars = len(data)
    bar_width = total_width / n_bars
    bars = []
    for i, (name, values) in enumerate(data.items()):
        x_offset = (i - n_bars / 2) * bar_width + bar_width / 2
        for x, y in enumerate(values):
            bar = ax.bar(x + x_offset, y, width=bar_width * single_width, color=sns.color_palette("deep")[i % len(colors)])
        bars.append(bar[0])
    if legend:
        ax.legend(bars, data.keys(), fontsize = args.legend_size)
    if legend == 2:
        ax.legend(bars, data.keys(), fontsize = args.legend_size,bbox_to_anchor=(1, 0.5))


def draw_hist(args):

    dict_mutant = get_data(args)
    if not args.hist_context:
        meth_TE = {"CG" : [m[0] for m in dict_mutant.values()],
        "CHG" : [m[1] for m in dict_mutant.values()],
        "CHH" : [m[2] for m in dict_mutant.values()]}
        fig, ax = plt.subplots(figsize = (16,10))
        ax.yaxis.offsetText.set_fontsize(20)
        bar_plot(ax, meth_TE, args, total_width=.8, single_width=.9)
        plt.xticks(np.arange(len(dict_mutant.keys())),[key.upper() for key in dict_mutant.keys()], fontsize = args.xsize)
        plt.yticks(fontsize = args.yaxis_size)
        plt.ylabel(args.ylab_hist, fontsize = args.ysize)
        plt.title(args.title_hist, fontsize = 25)
        plt.savefig(args.output_hist, bbox_inches='tight')
    else:
        width_mid = 0.8/len(dict_mutant) / 2
        for i,context in enumerate(["CG","CHG","CHH"]):
            context_TE = {name:[m[i]] for name, m in dict_mutant.items()}
            fig, ax = plt.subplots(figsize = (16,10))
            ax.yaxis.offsetText.set_fontsize(20)
            bar_plot(ax, context_TE, args, total_width=.8, single_width=.9, legend = False)
            plt.xticks(np.arange(-0.4+width_mid,0.8-width_mid,width_mid)[::2][:len(dict_mutant)],[key.upper() for key in dict_mutant.keys()], fontsize = args.xsize)
            plt.yticks(fontsize = args.yaxis_size)
            plt.ylabel(args.ylab_hist, fontsize = args.ysize)
            plt.title(re.sub("stretches",f"{context} stretches",args.title_hist), fontsize = 25)
            plt.savefig(re.sub("CG",context,args.output_hist), bbox_inches='tight')

def draw_violin(args):
    header = ["chr","start","end","read","strand","length", "pvalue_expo"]
    for context in ["CG","CHG","CHH"]:
        df = []
        with open(args.file) as filin:
            for line in filin:
                if line.startswith("#"):
                    continue
                items = line.split()
                if len(items) == 0:
                    continue
                name = items[0]
                path = items[1]
                context_path = re.sub("CG",context, path)
                stretch = pd.read_csv(context_path, names = header, sep="\t")
                stretch_length = stretch["end"] -  stretch["start"]
                stretch_length = np.log10(stretch_length)
                df.append(pd.concat([ pd.DataFrame(stretch_length, columns = ["length"]),
                                    pd.DataFrame([name]*stretch_length.shape[0], columns = ["mutant"])], axis = 1))
        data = pd.concat(df)

        fig, ax = plt.subplots(figsize = (16,10))
        sns.violinplot(x = "mutant",y = "length", data = data, palette = "deep")
        plt.xlabel("")
        plt.ylabel("Length in bp (log10)", fontsize = 20)
        plt.yticks(fontsize = args.yaxis_size)
        plt.xticks(fontsize = args.xsize)
        plt.title(re.sub("CG",context, args.title_violin), fontsize = 25)
        plt.savefig(re.sub("CG",context, args.output_violin),  bbox_inches='tight')
        plt.show()

def draw_boxplot(args):

    header_pvalue = ["chr","start","end","read","strand","length","pvalue_expo"]

    for context in ["CG","CHG","CHH"]:
        with open(args.file) as filin:
            for line in filin:
                if line.startswith("#"):
                    continue
                items = line.split()
                if len(items) == 0:
                    continue
                name = items[0]
                file = re.sub("CG",context,items[1])
                print(file)
                stretch_data = pd.read_csv(file, names = header_pvalue, sep="\t")
                print(stretch_data["length"].shape[0])
                if stretch_data["length"].shape[0] > 10000:
                    tick_spacing = 5
                else:
                    tick_spacing = 1

                plt.figure(figsize = (40,15))
                plt.xticks(fontsize = 30, rotation = 90)
                plt.yticks(fontsize = 30)
                plt.ylabel("-log10 pvalue", fontsize = 30, labelpad = 15)
                plt.xlabel("Stretch length in CMCs", fontsize = 30, labelpad = 15)
                plt.title(f"{name} {context} stretches", fontsize = 50)
                ax = sns.boxplot(x = stretch_data["length"],
                            y = stretch_data["pvalue_expo"])
                ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
                ax.hlines(-np.log10(0.05 / stretch_data.shape[0]), 0, stretch_data.length.unique().shape[0] - 1, color = "red", lw = 10)
                plt.savefig(f"{args.outdir}/{name}_{context}_stretch_boxplot.pdf", bbox_inches='tight')


if __name__ == "__main__":
    args = get_args()
    if args.hist:
        draw_hist(args)
    if args.violin:
        draw_violin(args)
    if args.box:
        draw_boxplot(args)