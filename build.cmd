@echo off

:main
set TARGET_NAME="meso-mpi"
set SRC_DIR=./
set CMAKE_DIR=./cmake-build-release

if not exist cmake-build-release mkdir cmake-build-release

cmake -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles" -S %SRC_DIR% -B %CMAKE_DIR%
cmake --build %CMAKE_DIR% --target clean
cmake --build %CMAKE_DIR% --target %TARGET_NAME% -- -j 8
goto:eof
