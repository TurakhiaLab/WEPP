# create build directory
startDir=$PWD
cd $(dirname "$0")
mkdir -p build
cd build
buildDir=$PWD

# install minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf -
cp minimap2-2.24_x64-linux/minimap2 .
chmod +x minimap2

# install samtools
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9; 
make
cd $buildDir
mv samtools-1.9/samtools samtools

# install bw2-mem2
#git clone https://github.com/bwa-mem2/bwa-mem2
#cd bwa-mem2
#git submodule init
#git submodule update
#git clone --recursive https://github.com/bwa-mem2/bwa-mem2
#cd bwa-mem2
#make

cd $startDir
