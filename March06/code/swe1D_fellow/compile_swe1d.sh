#   Shallow Water Equations 1D Simulator (SWE1D)
#  ===========================================================================================
#   Rectangular cross section
#  ===========================================================================================
#   Version 1.0 - Jan 2025
#  ===========================================================================================
#   Computational Hydraulics Group - University of Zaragoza   
#  =========================================================================================== 

#!/bin/bash

# Nombres de los archivos fuente
SRC1="swe1d.c"
SRC2="lib/shallow_water.c"
OUTPUT="runSWE1D"
resultspath="outputFiles"

rm $OUTPUT
rm -r $resultspath

# Compilar los archivos con g++
g++ -o $OUTPUT $SRC1 $SRC2 -lm -lboost_iostreams -lboost_system -lboost_filesystem
mkdir $resultspath

# Verificar si la compilación fue exitosa
if [ $? -eq 0 ]; then
    echo "Code compiled succesfully. Run with ./$OUTPUT"
else
    echo "Compilation fail"
fi
