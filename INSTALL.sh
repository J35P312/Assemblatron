git clone --recursive https://github.com/lh3/fermikit.git
cd fermikit
make
cd ..
cd htsbox
make
cd ..

cd TIDDIT
mkdir build
cd build
cmake ..
make
cd ..
cd ..

git clone https://github.com/lh3/seqtk.git
cd seqtk 
make
