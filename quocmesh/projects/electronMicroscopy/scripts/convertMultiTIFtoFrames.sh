#!/bin/sh
for FILE in *.tif
do
  BASE="`basename $FILE .tif`"
  convert ${FILE} TIF:frames/${BASE}%03d.tif
done
