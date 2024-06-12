@echo off

:main

cmake -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles" -S "./" -B "./cmake-build-windows"
cmake --build "./cmake-build-windows" --target clean
cmake --build "./cmake-build-windows" --target "meso-mpi" -- -j 8
goto:eof
