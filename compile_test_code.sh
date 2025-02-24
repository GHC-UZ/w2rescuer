#   Test code
#  ===========================================================================================
#   Version 1.0 - Jan 2025
#  ===========================================================================================
#   Computational Hydraulics Group - University of Zaragoza   
#  =========================================================================================== 

#!/bin/bash

# Nombres de los archivos fuente
SRC1="test_code.c"
OUTPUT="run_test_code"
resultspath="outputFiles"

rm $OUTPUT
rm -r $resultspath

# Compilar los archivos con g++
g++ -o $OUTPUT $SRC1 -lm -lboost_iostreams -lboost_system -lboost_filesystem
mkdir $resultspath

# Verificar si la compilación fue exitosa
if [ $? -eq 0 ]; then
    echo "Code compiled succesfully. Run with ./$OUTPUT"
else
    echo "Compilation fail"
fi
