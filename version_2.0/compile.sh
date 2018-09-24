#!/bin/bash

DEFAULT=/data1/home/stephan/Photoionization_SAE_PW/version_2.0
SOURCE=/data1/home/stephan/Photoionization_SAE_PW/version_2.0

if [[ -z $1 ]]
then
   DIRECTORY=$DEFAULT
else
   DIRECTORY=$1
fi


echo "$DIRECTORY"
icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -g -openmp  -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${MKLROOT}lib/intel64 -mkl -liomp5 -lpthread -ldl  ${DIRECTORY}/photoion_comp.cpp -o ${DIRECTORY}/photoionization_comp.exe
