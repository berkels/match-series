#!/bin/sh
export "QUOC_NO_SYSTEM_PAUSE=ON"
for FILE in *.tif
do
  ./convertTIFFToQuoc $FILE
done 
