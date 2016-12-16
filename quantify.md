# Quantify transcripts with RNA-seq

Salmon is a fast and accurate tool for quantifying trancript
abundances, which can then be used to also generate gene-level
quantifications.

A binary for running Salmon can be found here:

https://github.com/COMBINE-lab/salmon/releases

A workflow of how to run Salmon can be found here:

https://combine-lab.github.io/salmon/getting_started/

First you'll need to download or generate a FASTA file of all of the
transcript sequences. You can download a FASTA file of transcript
sequence from
[ENSEMBL](http://useast.ensembl.org/info/data/ftp/index.html) or
[GENCODE](https://www.gencodegenes.org/). You can also generate such a
file using a tool such as
[rsem-prepare-reference](http://deweylab.biostat.wisc.edu/rsem/rsem-prepare-reference.html)
with a GTF file and a genome FASTA file.

Then, the two Salmon steps are generating an index, and then quantifying each
sample against that index. Good to use a descriptive index name.

```
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa -i Homo_sapiens.GRCh38.cdna.all.fa_Salmon-0.7.2_index
```

(This takes ~6 minutes.)

Quantifying reads against this index might look like:

```
salmon quant -i Homo_sapiens.GRCh38.cdna.all.fa_Salmon-0.7.2_index \
  -l IU \
  -1 fastq/sample_1.fastq.gz \
  -2 fastq/sample_2.fastq.gz \
  -p 6 -o quants/sample
```

The `-l IU` parameter specifies we have paired-end reads facing inward
and an unstranded protocol. For other library types, see the Salmon website.

See the workflow above for an example of how to incorporate this into
a bash script across a number of samples.
