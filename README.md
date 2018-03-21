# Overview
This module lets one create one or more redundancy descriptors,
which then may be applied to a set of files distributed across a group of processes.

Currently, it assumes the parent group of processes is MPI_COMM_WORLD.

Usage is documented in src/redset.h.

It defines two interfaces: a low-level "redset" API and a higher level "ER" interface.

# Building

To build KVTree:

    git clone git@github.com:LLNL/KVTree.git KVTree.git

    mkdir build
    mkdir install
    
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=../install -DMPI=ON ../KVTree.git
    make clean
    make
    make install
    make test

To build redset:

    cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_KVTREE_PREFIX=`pwd`/install .

# Testing
Some simple test programs exist in the test directory.

To build a test for the redset API:

    mpicc -g -O0 -o test_redset test_redset.c -I../install/include -L../install/lib64 -lkvtree -I../src -L../src -ler

To build a test for the ER API:

    mpicc -g -O0 -o test_er test_er.c -I../install/include -L../install/lib64 -lkvtree -I../src -L../src -ler
