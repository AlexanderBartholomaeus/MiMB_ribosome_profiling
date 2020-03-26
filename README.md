# Codon resolution analysis of ribosome profiling data
Here present a workflow of precise positioning of ribosome profiling data.

## Introduction
Translation is a central biological process in living cells. The ribosome profiling approach [Ingolia et al.](https://www.ncbi.nlm.nih.gov/pubmed/19213877) revolutionized the field and enables assessing translation on a global, cell-wide level. Extracting versatile information from the ribosome profiling data requires specialized expertise for handling the sequencing data that is not available to the broad community. Here, we provide an easy-to-use and modifiable workflow that enables precise positioning of the ribosome-protected fragments or translating ribosomes for determining codon-specific translation features. The workflow is complemented with simple step-by-step explanations and is published in [Methods in Molecular Biology - link will come soon](). 

## Software and Installation

### Required software
The pipeline was testen on Linux (Ubuntu 18.04. LTS). The following software is required:

- [R](https://www.R-project.org)
- [bedtools](https://bedtools.readthedocs.io)
- [samtools](http://www.htslib.org/)
- optional for mapping:
  - [cutadapt](https://cutadapt.readthedocs.io)
  - [bowtie](http://bowtie-bio.sourceforge.net)
 
### Installation notes
There different ways of installing the different tools. For most of the tools we provide a native and installation via [bioconda](http://bioconda.github.io/index.html). Bioconda allows easy install of many tools. Find infos on how to set it up on the [conda install page](http://bioconda.github.io/user/install.html#install-conda).

#### R
R is available for most Linux systems and can be installed via:

```
# ubuntu install (root rights required)
sudo apt-get install R

```

After installing R the following to R packages are needed:

- foreach
- doParallel

You can run the `install.R` script or install it manually in R with following command:
```
install.packages(c('foreach','doParallel'))
``` 

#### bedtools 
Bedtools is available for most Linux systems and can be installed via:

```
# ubuntu install (root rights required)
sudo apt-get install bedtools

```

Also a conda installation is possible:

```
conda install bedtools
```

To install the binaries visit the [bedtools installation page](https://bedtools.readthedocs.io/en/latest/content/installation.html).

#### samtools
Samtools can be installed via bioconda:
```
conda install samtools
```

or must be build from source [samtools install page](http://www.htslib.org/download/)

#### cutadapt
Cutadapt can be installed via the Python manager:

```
# install as local user
python3 -m pip install --user --upgrade cutadapt

```

or via conda:

```
conda install cutadadpt
```

#### Bowtie
Bowtie can be downloaded from [bowtie webpage](http://bowtie-bio.sourceforge.net) as a binary or install via conda:

```
conda install bowtie
```

## Example data 

We use publically available data published in [Mohammad et al.](https://elifesciences.org/articles/42591) to illustrate the use of the workflow. For the example data file `rpf.fastq.gz` we use the first one million reads of [SRR1734437](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1734437] which contains RPFs of *E. coli* MG1655. The example data can be found in the `example_data` folder.

The example data includes also prepared genome and annotation data. For the optional mapping a reference is needed, we provide the genome here. For the other steps a BED file with CDS information of is needed. For eukaryotic organism different exons must have the same name (column 4).

## Example code using the example data to perform ribosome profiling

The full workflow with a simple step-by-step and comprehensive explanations and is published in [Methods in Molecular Biology - link will come soon](). Here we give the minimal version. 

To run the command some configuration files must be set. There are CSV files where you can simple modify the fields using a text editor or excel (save as CSV not as EXCEL). The files are in the `config` folder:

* cores_max.csv: this file set the maximum number of parallel cores to use
* mapping_config.csv: this file contains parameters to set for adapter cutting and mapping
* tools_locations.csv: this file contains the path/location of the used tools

The example commands expect that you execute the scripts from the terminal/console inside the cloned git folder. When you open a terminal/console you can navigate to the folder by the following command:

```
# change PATH by the path and folder structure of your machine
cd PATH/MiMB_ribosome_profiling

```

For the example commands `../out` is the output folder which is located one directory over the folder with the scripts. You can change the path and name if you like. The example code use RPF and (m)RNA read files. The (m)RNA files are used as a control because we expect a more uniform distribution of read length and no enrichment of read around the translation start and stop.

### [Optional] Sequencing read pre-processing and mapping

This step is optional, as you can used your own mapped data. It performs adapter trimming using `cutadapt` and mapping using `bowtie` (version 1 is stated to be more sensitive for short reads) using the `config/mapping_config.csv` to control parameters. 

```
# perform adapter cutting and mapping
Rscript mapping.R example_data/sequencing_data/ ../out example_data/genome_data/E_coli_genome.fa
``` 

### Statistics and filtering ribosomal reads by length

Here the read length distribution is generated to get an overview of the occurence of each read length.

```
# get read distribution info and plots
Rscript read_length_distribution.R ../out/mapping/ ../out
``` 

### RPF coverage within start and stop codon regions

The scripts used here select good expressed genes, prepare annotation files (to include 50 nucleotides upstream and downstream that start and stop) and generate a coverage plot using the full length of the reads.

```
# determine the genes with good expression levels
Rscript highest_expressed_genes.R ../out/mapping/rpf.sort.bam example_data/genome_data/E_coli_genes.bed ../out/ 0.1
# prepare the BED file
Rscript prepare_coverage.R ../out/highest_expressed_genes/highest_expressed_genes.bed
# generate plot to show coverage of RPFs
Rscript coverage_start_stop.R ../out/mapping/ ../out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed ../out/
```

### Calibration of RPFs

Here the precise calibration of the RPF reads is performed. After executing `Rscript calibration_count_plot.R ...` the coverage plots for each read length can be found in the calibration folder. For each read length there is a summary of reads in the certain frame (page 1) which should show a clear non-uniform distribution. Otherwise these read length should be discarded. The coverage plots (page 2) show the 5' assigned reads. You may change this to 3' (using `3prime` as a parameter). You have to manually determine the offset. Using the example data 24nt and 27nt reads look very good and the offset to calibrate for ribosomal P-site would be 9 and 12 nucleotides. This offsets has to be entered into the `calibration_config.csv` file in the `../out/calibration` folder. Subsequently, reads are calibrated, merged into one file and coverage plot using calibrated reads is performed.

```
# split reads by length, 5' assignment for read length of 24 to 30 nucleotides
Rscript split_by_length.R ../out/mapping/ ../out/ 5prime 24-30
# generate the plot to calibrate the RPFs
Rscript calibration_count_plot.R ../out/calibration/split_by_length/ ../out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed ../out/
# calibrate the reads accoring to manually determined offsets using the calibration_config.csv in the ../out/calibration folder
Rscript calibrate_reads.R ../out/mapping/ example_data/calibration_5prime_config.csv ../out
# merge different read length after calibration
Rscript merge_reads.R ../out/calibration/calibrated/
# generate a coverage plot using precisely calibrated reads
Rscript coverage_start_stop.R ../out/calibration/calibrated/ ../out/highest_expressed_genes/highest_expressed_genes_plus_50nt.bed ../out/calibrated_coverage
```
