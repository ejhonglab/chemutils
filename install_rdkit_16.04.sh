#!/usr/bin/env bash


# can use this to check boost version is sufficient
# dpkg -s libboost-dev | grep 'Version'

#sudo apt-get install flex bison build-essential python-numpy cmake python-dev sqlite3 libsqlite3-dev libboost-dev  libboost-python-dev libboost-regex-dev

cd ~/src
#git clone git://github.com/rdkit/rdkit.git

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

cmake -DRDK_BUILD_INCHI_SUPPORT=ON -DPYTHON_LIBRARY=/usr/lib/python3.5/config-3.5m-x86_64-linux-gnu/libpython3.5m-pic.a -DPYTHON_INCLUDE_DIR=/usr/include/python3.5/ -DPYTHON_EXECUTABLE=/usr/bin/python3 ..
make -j4
make install

ctest
