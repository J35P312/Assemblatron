git clone --recursive https://github.com/lh3/fermikit.git
cd fermikit
make clean
make
cd ..

pip install BESST

cd scripts
python setup.py build_ext --inplace

cd ..

git clone https://github.com/ablab/quast.git
cd quast
pip install -e .
cd ..

git clone https://github.com/lh3/bfc.git -b kmc
cd bfc
make ..

conda install -c bioconda kmc
