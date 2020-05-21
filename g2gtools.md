# Setting up g2gtools on longleaf cluster

```
module load anaconda
git clone https://github.com/churchill-lab/g2gtools.git
cd g2gtools
git checkout develop
conda create -n g2gtools python=2
conda activate g2gtools
conda install numpy pysam cython natsort bx-python future
python setup.py install
```
# Running g2gtools

Refer to <http://churchill-lab.github.io/g2gtools/#overview> for more detail.

First active the g2gtools environment 
```
module load anaconda
conda activate g2gtools
```
- Create VCI file from VCF (use --diploid for diploid analysis)
```
g2gtools vcf2vci -i <chr1.VCF> -f <chr1.fasta> -s <strain name> -o <output.vci> 
```
- Incorporate SNPs into the reference genome
```
g2gtools patch -i <chr1.fasta> -c <output.vci> -o <patch.fa>
```
- Incorporate indels into the reference genome and get custom diploid genome for the sample
```
g2gtools transform -i <patch.fa> -c <output.vci> -o <transform.fa>
```
- Liftover gene annotation onto sample coordinates
```
g2gtools convert -c <output.vci> -i <reference.coordinates.gtf> -o <sample.gtf>
```
- Parse custom gene annotation into a G2G database file 
```
g2gtools gtf2db -i <sample.gtf> -o <sample.db>
```
Note the process time is approximately 1 min/100,000 records. It is recommended to submit a bash file if the record number is large
The record number can be checked by using 
```
wc -l <sample.gtf>
```


