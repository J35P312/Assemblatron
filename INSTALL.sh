git clone --recursive https://github.com/lh3/fermikit.git
cd fermikit
make clean
make
cd ..

git clone https://github.com/J35P312/TIDDIT.git
cd TIDDIT
mkdir build
cd build
cmake ..
make
cd ..
cd ..

pip install BESST
git clone https://github.com/J35P312/SVDB.git
cd SVDB
pip install -e .
cd ..
