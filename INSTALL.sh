git clone --recursive https://github.com/lh3/fermikit.git
cd fermikit
make clean
make
cd ..

pip install BESST

cd scripts
python setup.py build_ext --inplace
