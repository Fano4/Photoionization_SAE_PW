#!/bin/bash

DEFAULT=${HOME}/Photoionization_SAE_PW/version_1.0
SOURCE=${HOME}/Photoionization_SAE_PW/version_1.0

if [[ -z $1 ]]
then
   DIRECTORY=$DEFAULT
else
   DIRECTORY=$1
fi


echo "$DIRECTORY"
icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -openmp  -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${MKLROOT}lib/intel64 -mkl -liomp5 -lpthread -ldl  ${SOURCE}/photoion_comp.cpp -o ${DIRECTORY}/PICE_2PI.exe
