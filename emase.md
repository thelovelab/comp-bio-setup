# Setting up EMASE on longleaf cluster

```
module load anaconda
conda create --name emase
conda activate emase 
conda install -c anaconda numpy==1.8.2 scipy==0.13.3 pysam pytables biopython cython 
conda install -c kbchoi emase
git clone git@github.com:churchill-lab/emase.git
cd emase
git checkout release/0.10.17
python setup.py install
./scripts/run-emase
```
