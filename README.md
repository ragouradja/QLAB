# QLAB

## Conda environment

Clone this repository : 

    $ git clone git@github.com:ragouradja/QLAB.git
Or : 

    $ git clone https://github.com/ragouradja/QLAB.git
Then :

    $ cd QLAB

Create the conda environment : 
    
    $ conda env create -f env/env.yml
Load the conda environment : 
    
    $ conda activate ont_env

You are ready to run any scripts. Some test data are available in `/mnt/data5/rradjas/test_script`.

## Knowns errors

If all fast5 failed during Tombo or DSP, you need to add `export HDF5_PLUGIN_PATH...` variable in your bashrc with the correct path to hdf5 plugin that you will download :

```bash
wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
tar zxvf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
export HDF5_PLUGIN_PATH=/abslolute/path/to/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
```

In some scripts, the `sort` can fails due to large file. You must add the option `sort -T tmp bigfile` but sometimes it still doesn't work because he is not finding the `tmp` folder :

`sort: cannot create temporary file in 'tmp': No such file or directory`

To avoid this, the tmp folder is created before the sort command in the scripts. You can remove any tmp folder that have been created with scripts.

To compress large file faster than classical compression method, I use `pigz` which will use multiprocessing to compress file :
https://rachaellappan.github.io/pigz/

which can be installed with
```bash
sudo apt install pigz 
```
