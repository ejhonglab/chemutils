#!/usr/bin/env bash


# can use this to check boost version is sufficient
# dpkg -s libboost-dev | grep 'Version'

#sudo apt-get install flex bison build-essential python-numpy cmake python-dev sqlite3 libsqlite3-dev libboost-dev  libboost-python-dev libboost-regex-dev

cd ~/src
git clone git://github.com/rdkit/rdkit.git

printf "\n\n# For rdkit cheminformatics library\n" >> ~/.bashrc
printf "export RDBASE=\$HOME/src/rdkit\n" >> ~/.bashrc
printf "export LD_LIBRARY_PATH=\$RDBASE/lib:\$LD_LIBRARY_PATH\n" >> ~/.bashrc
printf "export PYTHONPATH=\$RDBASE:\$PYTHONPATH\n\n" >> ~/.bashrc
. ~/.bashrc

cd rdkit/External/INCHI-API/
./download-inchi.sh
cd ../..

mkdir build
cd build

# TODO note: before cmake, may need to delete coordgen[libs] + coordgen*.tar.gz
# + maeparser under External/CoordGen, as in:
# https://github.com/rdkit/rdkit/issues/2432

# (could never get it to work w/ anaconda)
# For one install of Anaconda, these paths seemed right:
# -DPYTHON_LIBRARY=/home/tom/anaconda3/lib/python3.7/config-3.7m-x86_64-linux-gnu/libpython3.7m.a
# -DPYTHON_INCLUDE_DIR=/home/tom/anaconda3/include/python3.7m/
# -DPYTHON_EXECUTABLE=/home/tom/anaconda3/bin/python3
# Some version of this flag may be useful to fix inability to find boost
# headers:
# -DCMAKE_CXX_FLAGS=-isystem\ /usr/include 

cmake -DRDK_BUILD_INCHI_SUPPORT=ON -DPy_ENABLE_SHARED=1 -DPYTHON_LIBRARY=/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5m-pic.a -DPYTHON_INCLUDE_DIR=/usr/include/python3.5/ -DPYTHON_EXECUTABLE=/usr/bin/python3 ..
make -j4
make install

ctest
