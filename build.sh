#!/bin/bash

cmake -DCMAKE_BUILD_TYPE=Relaese -G "Unix Makefiles" -S "./" -B "./cmake-build-linux"
cmake --build "./cmake-build-linux" --target clean
cmake --build "./cmake-build-linux" --target "meso-mpi" -- -j 18
