#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Build and install TorchFort
cd /tmp/sachinbm/TorchFort/build/
cmake -DCMAKE_INSTALL_PREFIX=/tmp/sachinbm/TorchFort/build/install \
      -DCMAKE_C_COMPILER=$CC \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_Fortran_COMPILER=$FC \
      -DTORCHFORT_BUILD_EXAMPLES=1 \
      -DTORCHFORT_BUILD_TESTS=1 \
      -DTORCHFORT_ENABLE_GPU=0 \
      -DTORCHFORT_YAML_CPP_ROOT=$TORCHFORT_YAML_CPP_ROOT \
      -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH;$(python -c 'import torch; print(torch.utils.cmake_prefix_path)')" \
      ..
make -j1 install VERBOSE=1

# Configure and build Neko with TorchFort
cd /tmp/sachinbm/neko-tf/neko/
./regen.sh
env FCFLAGS="-I/tmp/sachinbm/TorchFort/build/install/include" \
    LDFLAGS="-L/tmp/sachinbm/TorchFort/build/install/lib64" \
    LIBS="-ltorchfort_fort -ltorchfort" \
    ./configure FC=$FC CC=$CC MPIFC=$MPIFORT MPICC=$MPICC --prefix="/tmp/sachinbm/neko-tf/neko/install"
make -j4 install

# Build and run the example
cd examples/turb_channel
/tmp/sachinbm/neko-tf/neko/install/bin/makeneko turb_channel.f90
./neko les_2.case
