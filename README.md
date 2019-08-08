# Overview

[![Build Status](https://api.travis-ci.org/ECP-VeloC/redset.png?branch=master)](https://travis-ci.org/ECP-VeloC/redset)

This module lets one create one or more redundancy descriptors,
which then may be applied to a set of files distributed across a group of processes.

Usage is documented in [src/redset.h](src/redset.h) and [doc/rst/redset.rst](doc/rst/redset.rst).
Also see the example program in the [test](test) directory.

Implementation details can be found in [doc/rst/implementation.rst](doc/rst/implementation.rst) and [doc/rst/schemes.rst](doc/rst/schemes.rst).

# Building

To build dependencies:

    git clone git@github.com:ECP-VeloC/KVTree.git  KVTree.git
    git clone git@github.com:ECP-VeloC/rankstr.git rankstr.git

    mkdir install

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=../install -DMPI=ON ../KVTree.git
    make clean
    make
    make install
    make test
    cd ..

    rm -rf build
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=../install -DMPI=ON ../rankstr.git
    make clean
    make
    make install
    cd ..

To build redset:

    cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=./install -DWITH_KVTREE_PREFIX=`pwd`/install -DWITH_RANKSTR_PREFIX=`pwd`/install .
    make
    make install

# Testing
Some simple test programs exist in the test directory.

To build a test for the redset API:

    mpicc -g -O0 -o test_redset test_redset.c -I../install/include -L../install/lib64 -lkvtree -lrankstr -I../src -L../src -lredset

## Release

Copyright (c) 2018, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
<br>
Copyright (c) 2018, UChicago Argonne LLC, operator of Argonne National Laboratory.


For release details and restrictions, please read the [LICENSE]() and [NOTICE]() files.

`LLNL-CODE-751725` `OCEC-18-060`
