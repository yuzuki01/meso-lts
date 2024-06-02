#!/bin/bash

make() {
  SRC_DIR=./src/$TARGET_NAME
      CMAKE_DIR=./cmake-build-release/$TARGET_NAME

      if [ ! -d "cmake-build-release" ]; then
          mkdir cmake-build-release
      fi
      if [ ! -d "$CMAKE_DIR" ]; then
          mkdir "$CMAKE_DIR"
      fi

      cmake -DCMAKE_BUILD_TYPE=Debug -G "Unix Makefiles" -S "$SRC_DIR" -B "$CMAKE_DIR"
      cmake --build "$CMAKE_DIR" --target clean
      cmake --build "$CMAKE_DIR" --target "$TARGET_NAME" -- -j 8
}

case "$1" in
    core)
        TARGET_NAME=core
        make
        ;;
    mesh)
        TARGET_NAME=mesh
        make
        ;;
    solver)
        TARGET_NAME=solver
        make
        ;;
    *)
        TARGET_NAME=core
        make
        TARGET_NAME=mesh
        make
        TARGET_NAME=solver
        make
        TARGET_NAME=meso
        make
        ;;
esac
