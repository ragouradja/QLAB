# Histogram, Violin and Box plot for stretches data

```bash
$ python stretch_plot.py --help

usage: stretch_plot.py [-h] --file FILE [--outdir OUTDIR] [--hist] [--box] [--violin] [--hist_context] [--title_hist]
                       [--title_violin] [--ylab_hist] [--ylab_violin] [--ylab_box] [--output_hist] [--output_violin]
                       [--xsize] [--ysize] [--yaxis_size] [--legend_size]

Plotting histogram

optional arguments:
  -h, --help            show this help message and exit
  --file FILE, -f FILE  Path to file containing path for CG stretch file for each mutant (1 mutant per line : NAME
                        PATH) (default: None)
  --outdir OUTDIR, -o OUTDIR
                        Path to output folder (default: images)
  --hist, -hi           Draw histogram (default: False)
  --box, -b             Draw boxplot (default: False)
  --violin, -vi         Draw violin plot (default: False)
  --hist_context, -hc   Plot histogram per context (default: False)
  --title_hist , -th    Title of the histogram (default: Normalized number of significant stretches in ONT mutants)
  --title_violin , -tv
                        Title of the violin plot (default: Length of statistically significant CG stretches)
  --ylab_hist , -yh     Label of y axis (default: Normalized amount)
  --ylab_violin , -yv   Label of y axis (default: Length in bp (log10))
  --ylab_box , -yb      Label of y axis (default: -log10 pvalue)
  --output_hist , -outh
                        Output filename of histogram with extension (pdf, png...) (default: hist_stretch.pdf)
  --output_violin , -outv
                        Output filename of violin plot with extension (pdf, png...) (default: CG_stretch_violin.pdf)
  --xsize , -xs         Size of xaxis label (default: 15)
  --ysize , -ys         Size of yaxis label (default: 15)
  --yaxis_size , -yts   Size of yaxis values (default: 15)
  --legend_size , -ls   Size of legend (default: 15)
```

The `--file` argument is mandatory and need to be a file containing : `NAME | PATH | NORM`
* `NAME` is the name that should appear in the plot for this sample
* `PATH` is the path to CG stretch data. CHG and CHH files will be automatically found if they are in the same folder as CG file.
* `NORM` contains a value (average coverage for each sample for example) which will be used to divide the total number of stretches (and normalize them). This column is required only if you use `--hist`, otherwise it is not considered.

`NAME` and `PATH` are mandatory, `NORM` can be empty. Example with `NORM` value :
 
```bash
$ cat data_path.txt
WT test_data/WT/CG_genome_stretch.bed 53
CMT3 test_data/CMT3/CG_genome_stretch.bed 16.23
CMT2 test_data/CMT2/CG_genome_stretch.bed 11.8
RDR2 test_data/RDR2/CG_genome_stretch.bed 15.6
DRM12CMT3 test_data/DRM12CMT3/CG_genome_stretch.bed 29.57
```


## Histogram and Violin plot
To use `--hist` `--violin`, it is preferable to use the files containing statistically significant stretches files : 

```bash
$ python stretch_plot.py --file path_stretch_sign.txt --hist --violin
```

If you want the histogram to be per context, you can use `--hist_context` argument : 

```bash
$ python stretch_plot.py --file path_stretch_sign.txt --hist --hist_context
New name of histogram output file : CG_hist.pdf
You can change it with --output_hist
Drawing histogram...
```
The default output filename for histogram is `hist_stretch.pdf` but if you ask for one histogram per context, `CG` must be on the output filename in order to create a pdf file for `CG` `CHG` and `CHH`, otherwise a default name is used with `CG` in the name. You can change the name of the output file with `--output_hist`.

## Boxplot of pvalues
To use `--box`,  it is preferable to use the files containing all stretches files to see stretch that are not statistically significant and those that are : 

```bash
$ python stretch_plot.py --file path_stretch.txt --box
```

## Title, axis name...
You can change the title, x_axis or y_axis name (and fontsize) with other argument available. For example, change the output filename, title, ylab and ylab size :

```bash
$ python stretch_plot.py --file path_stretch_sign.txt --hist --title_hist "New title" --ylab_hist "New Y lab" --ysize 50 --output_hist new_hist.pdf
```


