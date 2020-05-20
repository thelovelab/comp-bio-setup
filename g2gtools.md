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
