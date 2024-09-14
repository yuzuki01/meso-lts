#!/bin/bash

DATE=$(date +"%Y%m%d")

ZIPFILE="meso-$DATE.zip"

FILES_TO_ZIP="./src ./files ./include ./lib ./build.sh ./CMakeFiles.txt ./README.md"

zip -r "$ZIPFILE" $FILES_TO_ZIP

echo "Zip files to: $ZIPFILE"
