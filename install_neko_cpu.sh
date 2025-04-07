#!/bin/bash
set -e

# Define path variables
TORCHFORT_INSTALL="/tmp/sachinbm/TorchFort/build/install"
OPENMPI_INSTALL="/tmp/sachinbm/openmpi-5.0.5/build/install"
NEKO_PREFIX="/tmp/sachinbm/neko-tf/neko/install"
GFORTRAN="/usr/bin/gfortran-13"
GCC="/usr/bin/gcc-13"
MPIFORT="${OPENMPI_INSTALL}/bin/mpifort"
MPICC="${OPENMPI_INSTALL}/bin/mpicc"

# Navigate to the neko directory and regenerate files
cd /tmp/sachinbm/neko-tf/neko
./regen.sh

# Configure using environment variables with the specified paths
env FCFLAGS="-I${TORCHFORT_INSTALL}/include" \
    LDFLAGS="-L${TORCHFORT_INSTALL}/lib64" \
    LIBS="-ltorchfort" \
    ./configure FC=${GFORTRAN} CC=${GCC} \
    MPIFC=${MPIFORT} MPICC=${MPICC} \
    --prefix=${NEKO_PREFIX}

# Build and install
make install
