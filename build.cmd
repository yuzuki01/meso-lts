@echo off

:main
if "%1" == "core" (
    set TARGET_NAME=core
    call:make
    goto:eof
)
if "%1" == "mesh" (
    set TARGET_NAME=mesh
    call:make
    goto:eof
)
if "%1" == "solver" (
    set TARGET_NAME=solver
    call:make
    goto:eof
)


set TARGET_NAME=core
call:make
set TARGET_NAME=mesh
call:make
set TARGET_NAME=solver
call:make
set TARGET_NAME=meso
call:make
goto:eof

:make
set SRC_DIR=./src/%TARGET_NAME%
set CMAKE_DIR=./cmake-build-release/%TARGET_NAME%

if not exist cmake-build-release mkdir cmake-build-release
if not exist cmake-build-release/%TARGET_NAME% mkdir cmake-build-release/%TARGET_NAME%

cmake -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles" -S %SRC_DIR% -B %CMAKE_DIR%
cmake --build %CMAKE_DIR% --target clean
cmake --build %CMAKE_DIR% --target %TARGET_NAME% -- -j 8
goto:eof
