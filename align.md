We can align reads to the genome using a software such as 
[hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml).

Download the appropriate index from the hisat2 website to our 
annotation directory, and then run `tar -xvf` on the file that 
you download (it will create a directory).

A simple bash script for running hisat2 on some paired-end reads
is then:

```
#!/bin/bash
#
#SBATCH --job-name=align
#SBATCH --time=120
#SBATCH --mem=10000
#SBATCH -n 12
#SBATCH -N 1

module load hisat2
module load samtools

hisat2 -p 12 -x hisat2_grch37_snp/genome_snp \
  -1 reads_01.fastq.gz \
  -2 reads_02.fastq.gz \
  -S file.sam

samtools view -bS file.sam > file.bam
samtools sort file.bam -o file.sorted.bam
samtools index file.sorted.bam
```

Note: if you are generating simulated reads with Polyester, 
you will need to add `-f` for FASTA input.
