#!/bin/bash

base1="${1::${#1}-8}"
base2="${1::${#1}-14}"

echo ${base1}
echo ${base2}

gdal_calc.py --calc="logical_and(logical_and(A==255, B==0), C==0)" -A $1 --A_band=1 -B $1 --B_band=2 -C $1 --C_band=3 --outfile=${base1}.mask.tif

gdal_translate -ot UInt16 -a_nodata 256 ${base2}RED-BLUE.rgb.tif ${base1}.acd.tif -co COMPRESS=LZW -b 1 -b 3 -b 3
