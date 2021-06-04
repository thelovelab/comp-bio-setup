# FASTQC and MultiQC

FASTQC and MultiQC are two tools for performing quality control (QC)
of sequence read data (fastq files). The most important things to
know:

* FASTQC is a program that is run on individual `fastq` or `fastqc.gz`
  files. It compiles basic statistics about read quality and base pair
  frequencies. Note that some of the quality checks are more
  appropriate for DNA-seq than RNA-seq. So just because you see a
  warning from a FASTQC report doesn't mean there is a problem.
* MultiQC is a program that aggregates statistics across many files,
  and many programs. If you run this program in a directory, it will
  search for files it recognizes and build an HTML report that
  summarizes all the information.
  
Both programs are available on the cluster. Just run:

```
module load fastqc
module load multiqc
```

How to incorporate into an RNA-seq workflow, see the following
Snakemake file, which will perform FASTQC, Salmon quantification, and
MultiQC report building all together:

<https://gist.github.com/mikelove/5a8134e57f652f970f1a176efc900cbe#file-snakefile_with_qc>

You could then run this program from the command line with, e.g.:

```
snakemake -j 4 --latency-wait 30 --cluster "sbatch -n 12 -N 1 --mem=10000 --time 60"
```

which would submit 4 jobs at one time to the cluster, each with 12 cores.
