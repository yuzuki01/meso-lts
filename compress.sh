#!/bin/bash

DATE=$(date +"%Y%m%d")

ZIPFILE="meso-$DATE.zip"

mkdir -p tmp/meso-mpi
cp -r ./src ./files ./include ./lib ./build.sh ./CMakeLists.txt ./README.md tmp/meso-mpi/
cd tmp || exit
zip -r "../$ZIPFILE" meso-mpi
cd ..
rm -rf tmp

echo "Zip files to: $ZIPFILE"
