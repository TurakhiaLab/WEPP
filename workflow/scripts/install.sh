rm -rf build
mkdir build
cd build
git clone https://github.com/microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install jsoncpp
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_DIR=${PWD}/oneTBB-2019_U9 -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake ..
# cmake -DCMAKE_BUILD_TYPE=Debug -DTBB_DIR=${PWD}/oneTBB-2019_U9 -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake -DCMAKE_TOOLCHAIN_FILE=${PWD}/vcpkg/scripts/buildsystems/vcpkg.cmake ..
