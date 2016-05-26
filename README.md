# CNA Processing

Purpose 
---
Process CNA data into formats usable by MAGI, HotNet, and CoMEt. 

Requires
---
* Python 2.7 (tested with 2.7.9)
* Matplotlib (tested with 1.4.2)
* Numpy (tested with 1.8.2)
* Gene target list (see below for format specification)
* Gene location list (see below for format specification)
* GISTIC2.0 output, unzipped to a flat directory or in a tar file (directory structure does not matter for tar file). The following files are the required output files:
    * Amplified peak file (typically `table_amp.conf_99.txt`)
    * Deleted peak file (typically `table_del.conf_99.txt`)
    * Focal matrix file (typically `focal_data_by_genes.txt`)
    * Focal segment file (typically `focal_input.seg.txt`).
* OR alternate format GISTIC2.0 output, unzipped to a flat directory or in a tar file. The following files are the required files, note that MAGI output is not possible with the alternate GISTIC2.0 format as it is missing the Focal segment file. Call the script with the `-ca` argument to check automatically for the alternate format:
    * Amplified peak file (typically `amp_genes.conf_99.txt`)
    * Deleted peak file (typically `del_genes.conf_99.txt`)
    * Focal matrix file (typically `focal_data_by_genes.txt`)

Setup
---
None required, beyond having necessary python libraries and gene target/location lists.

Quick Start
---
This will work assuming the GISTIC2.0 output data is using the standard naming 
convention. 

`python gistic2processing.py -g <path to gene dictionary> -d <path to data folder/tar file> -tg <path to gene target file>`

Output
---

#### HotNet2
A tab separated file, each row is a sample and associated gene mutations. The first column is the sample name, each entry after is a gene with either (A) for amplified or (D) for deleted appended on the end. The name is {prefix}_hotnet2.tsv. The prefix is a optional argument. If no prefix is supplied, the default is 'output'.

#### CoMEt
The same as HotNet2, named {prefix}\_all\_cna\_comet.tsv. There is an additional file named {prefix}\_name\_map\_comet.tsv, which provides a mapping from a shortened, more human-readable concatenated gene list to the full list of genes found in each peak. That file is a series of rows with the shortened name in the first column and the longer name in the second, tab separated.

#### MAGI
A five column, tab separated file with a header. Each row (in order) consists of a gene, a sample, the CNA type (Del/Amp), left position, and right position.


Usage
---
All options can be set by command line or via a configuration file. The 
command line will override configuration file settings, and the names and usage for
command line and configuration file options are the same. 

### Required ###

The following options must be specified by the user either at the command line or in the config file, and cannot be left to defaults:

Config/long argument | Short argument   | Input type | Description 
:-------------------------------------| :----- |:----- |:-----
gene_dictionary  | -g | Path to file |Path to the gene dictionary, a json file with a key/value pair of { gene : [[start location, end location, chromosome number]] }. See example folder in this repository.
data     | -d | Path to folder or tar file| Path to either a folder with a flat structure containing all of the GISTIC2.0 data, or the path to a tar file (does not have to have a flat internal file structure) with the GISTIC2.0 data.
target_genes      |  -tg | Path to file |Path to a file containing the target genes. File format is two columns, tab separated, with the gene name in the first column and Amp/Del/Both (not case sensitive) in the second column.
