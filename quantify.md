# Quantify transcripts with RNA-seq

Salmon is a fast and accurate tool for quantifying trancript
abundances, which can then be used to also generate gene-level
quantifications.

The latest release for running Salmon can be found here:

<https://github.com/COMBINE-lab/salmon/releases>

The linux binaries can be downloaded and run directly on the cluster. It's a good idea to just use the *latest* version of Salmon at the beginning of a project. Note that you have to create a new *index* if you use a new version of Salmon.

A workflow of how to run Salmon can be found here:

<https://combine-lab.github.io/salmon/getting_started/>

## Indexing the reference transcripts

First you'll need to download or generate a FASTA file of all of the
reference transcript sequences. You can download a FASTA file of transcript
sequence from
[GENCODE](https://www.gencodegenes.org/). You can also generate such a
file using a tool such as
[rsem-prepare-reference](http://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html)
with a GTF file and a genome FASTA file. It's a good idea to download
the GTF file whenever you download a transcripts FASTA file.

Then, the two Salmon steps are 1) generating an **index** (only needs to happen once), 
and then 2) **quant**ifying each sample against the index. Good to use a descriptive index name.
You should put the version number on the end of the index name (whatever 
version you are using, not 0.7.2).

Salmon indexing might look like this (this needs to be run either as a 
cluster job, see below, or from an interactive session).

```
salmon index -t gencode.vXX.transcripts.fa -i gencode.vXX_salmon_x.y.z
```

## Quantifying reads 

Quantifying reads with this index might look like:

```
salmon quant -i gencode.vXX_salmon_x.y.z \
  -l A --gcBias \
  -1 fastq/sample1_1.fastq.gz \
  -2 fastq/sample1_2.fastq.gz \
  -p 8 -o quants/sample1
```

Except, when working on the cluster we have to submit jobs to the cluster.
Therefore we write a bash script to submit to the cluster. However,
with bash scripts, we have to run a job for each sample, and it's
possible to make a mistake, e.g. if we change the sample name in the reads
but we forget to change the sample name of the output directory. 
Therefore, it's safer and a better idea to run many jobs *programmatically*:

I recommend using Snakemake to run Salmon across multiple files. My
Snakemake file for running Salmon is here:

<https://gist.github.com/mikelove/5a8134e57f652f970f1a176efc900cbe>

To run Snakemake on the cluster I use the following bash script
`snake.sh`:

```
#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=1440
#SBATCH --mem=1000

module load python/3.6.6
snakemake -j 4 --latency-wait 30 --cluster "sbatch --mem=10000 -N 1 -n 8"
```

The `-j 4` indicates to run 4 jobs at a time, each one with 10 Gb, and
with 8 threads for each job each.
