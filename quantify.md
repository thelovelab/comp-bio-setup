# Quantify transcripts with RNA-seq

Salmon is a fast and accurate tool for quantifying trancript
abundances, which can then be used to also generate gene-level
quantifications.

The latest release for running Salmon can be found here:

<https://github.com/COMBINE-lab/salmon/releases>

The Linux binaries can be downloaded and run directly on the
cluster. It's a good idea to just use the *latest* version of Salmon
at the beginning of a project. Note that you have to create a new
*index* if you use a new version of Salmon.

**Recommendation**: at the bottom of this file, I provide a
`Snakefile` which is the easiest way to run Salmon on the cluster
(everything is taken care of for you). So read over the following for
understanding what the steps look like, but for a project with many
files you shouldn't be running the steps one-by-one or with bash
scripts, but using Snakemake.

Also, a workflow from the Salmon authors on how to run Salmon can be
found here: 

<https://combine-lab.github.io/salmon/getting_started/>

## Indexing the reference transcripts

First you'll need to download or generate a FASTA file of all of the
reference transcript sequences. You can download a FASTA file of transcript
sequence from
[GENCODE](https://www.gencodegenes.org/). For organisms other than human or
mouse, you can use Ensembl FASTA 
(we recommend combining coding and non-coding RNA FASTA though).

(If GENCODE or Ensembl are not an option, you can also generate such a
file using a tool such as
[rsem-prepare-reference](http://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html)
with a GTF file and a genome FASTA file.)

It's a good idea to download the GTF file whenever you download a transcripts FASTA file.

Then, the two Salmon steps are 1) generating an **index** (only needs to happen once), 
and then 2) **quant**ifying each sample against the index. Good to use a descriptive index name.
You should put the Salmon version number on the end of the index name. The Salmon index will work
for that FASTA and any sample, as long as you use the same version of Salmon for quantification.
Re-indexing is no big deal, so when you upgrade Salmon, you need to re-index the transcripts FASTA.

Salmon indexing might look like this (this needs to be run either as a 
cluster job [best option], see below, or at least from an interactive session).

```
salmon index --gencode -t gencode.vXX.transcripts.fa -i gencode.vXX_salmon_x.y.z
```

(The `--gencode` is a special argument for GENCODE, which helps clean up the transcript names.)

## Quantifying reads 

Quantifying reads with this index might* look like this below
(note that the `quants` dir needs to be created first):

```
salmon quant \
  -l A --gcBias -p 12 \
  --numGibbsSamples 20 --thinningFactor 100 \
  -i gencode.vXX_salmon_x.y.z \
  -o quants/sample1 \
  -1 fastq/sample1_1.fastq.gz \
  -2 fastq/sample1_2.fastq.gz
```

* Except, when working on the cluster we want to submit jobs to the cluster!
Never run Salmon index or quant on the login node, and you don't want
to be sitting by your laptop during an interactive job quantifying 50 samples.
Therefore we could write a *bash script* to submit to the cluster. First I'll show 
what that would look like, then I'll give a better solution called Snakemake.

The batch script might look like:

```
#!/bin/bash
#
#SBATCH --job-name=quant
#SBATCH --time=60
#SBATCH --mem=10000
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mail-user=you@email.com
#SBATCH --mail-type=ALL

/path/to/salmon_x.y.z/bin/salmon quant \
  -l A --gcBias -p 12 \
  --numGibbsSamples 20 --thinningFactor 100 \
  -i gencode.vXX_salmon_x.y.z \
  -o quants/sample1 \
  -1 fastq/sample1_1.fastq.gz \
  -2 fastq/sample1_2.fastq.gz
```

And then you can submit this to cluster with `sbatch quant.sh`.

However, with bash scripts, we have to run a job for *each* sample, and it's
possible to make a mistake, e.g. if we change the sample name in the reads
but we forget to change the sample name of the output directory. 
Therefore, it's safer and a better idea to run many jobs *programmatically*
with [Snakemake](https://snakemake.readthedocs.io/en/stable/).

I recommend using Snakemake to run Salmon across multiple files. 

I have two Snakemake files for running Salmon, that can be found at
the link below.

1. `Snakefile` gives an example of how to run Salmon indexing and
   quantification over a number of files.
2. `Snakefile_with_QC` gives an example of running both Salmon
   indexing/quantification as well as QC with 
   [FASTQC and MultiQC](fastq_multiqc.md), all compiled into a single
   report. Note that to run this, you need to rename the file to
   `Snakefile` and load the FASTQC and MultiQC modules on the cluster.

<https://gist.github.com/mikelove/5a8134e57f652f970f1a176efc900cbe>

To run Snakemake on the cluster I use the following bash script
`snake.sh`:

```
#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=1440
#SBATCH --mem=1000

module load python
snakemake -j 4 --latency-wait 30 --cluster "sbatch --mem=10000 -N 1 -n 12"
```

The `-j 4` indicates to run 4 jobs at a time, each one with 10 Gb, and
with 12 threads for each job each.

# Index job

For completeness, a bash script for indexing:

```
#!/bin/bash
#
#SBATCH --job-name=index
#SBATCH --time=60
#SBATCH --mem=10000
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mail-user=you@email.com
#SBATCH --mail-type=ALL

/path/to/salmon_x.y.z/bin/salmon index -p 12 --gencode -t gencode.vXX.transcripts.fa -i gencode.vXX_salmon_x.y.z
```
