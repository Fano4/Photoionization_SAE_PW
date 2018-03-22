#!/bin/bash

DEFAULT=/CECI/home/ulg/cpt/svdwild/
SOURCE=/CECI/home/ulg/cpt/svdwild/photoionization_v1.0

if [[ -z $1 ]]
then
   DIRECTORY=$DEFAULT
else
   DIRECTORY=$1
fi


echo "$DIRECTORY"
icpp=icpc
${icpp} -openmp  -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${MKLROOT}lib/intel64 -mkl -liomp5 -lpthread -ldl  ${DIRECTORY}/photoionization_v1.0/photoion_comp.cpp -o ${DIRECTORY}/photoionization_comp.exe
