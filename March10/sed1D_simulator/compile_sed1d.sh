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
SRC1="sed1d.c"
SRC2="lib/shallow_water.c"
SRC3="lib/sediment.c"
OUTPUT="runSED1D"
resultspath="outputFiles"

##########################################################################################################################
DEBUG=no
#	valgrind --tool=memcheck --leak-check=full --track-origins=yes --show-reachable=yes
##########################################################################################################################


rm $OUTPUT
rm -r $resultspath

# Compilar los archivos con g++
if [ $DEBUG=yes ]; then
g++ -o $OUTPUT $SRC1 $SRC2 $SRC3 -g -lm -lboost_iostreams -lboost_system -lboost_filesystem
else
g++ -o $OUTPUT $SRC1 $SRC2 $SRC3 -lm -lboost_iostreams -lboost_system -lboost_filesystem
fi

# Verificar si la compilación fue exitosa
if [ $? -eq 0 ]; then
    echo "Code compiled succesfully. Run with ./$OUTPUT"
else
    echo "Compilation fail"
fi

mkdir $resultspath