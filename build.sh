#!/bin/bash

SRC_DIR="./"
CMAKE_DIR="./cmake-build-release"
TARGET_NAME="meso-mpi"

cmake -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles" -S "$SRC_DIR" -B "$CMAKE_DIR"
cmake --build "$CMAKE_DIR" --target clean
cmake --build "$CMAKE_DIR" --target "$TARGET_NAME" -- -j 8
