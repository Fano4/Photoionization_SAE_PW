#!/bin/bash

DEFAULT=/data1/home/stephan/Photoionization_SAE_PW/version_3.0
SOURCE=/data1/home/stephan/Photoionization_SAE_PW/version_3.0

GSL_ROOT=data1/home/stephan/gsl/

if [[ -z $1 ]]
then
   DIRECTORY=$DEFAULT
else
   DIRECTORY=$1
fi


echo "$DIRECTORY"
icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -g -openmp  -I/${GSL_ROOT}/include -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${GSL_ROOT}/lib -L/${MKLROOT}lib/intel64 -mkl -lgsl -liomp5 -lpthread -ldl -DMKL_Complex16="std::complex<double>"  ${DIRECTORY}/photoion_comp.cpp -o ${DIRECTORY}/photoionization_comp.exe
