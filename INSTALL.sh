git clone --recursive https://github.com/lh3/fermikit.git
cd fermikit
make
cd ..
cd htsbox
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
