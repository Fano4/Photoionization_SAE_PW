#!/bin/bash

DEFAULT=/data1/home/stephan/Photoionization_SAE_PW/version_3.0
SOURCE=/data1/home/stephan/Photoionization_SAE_PW/version_3.0

GSL_ROOT=data1/home/stephan/gsl
HF5_ROOT=data1/home/stephan/libraries/cpp/hdf5

if [[ -z $1 ]]
then
   DIRECTORY=$DEFAULT
else
   DIRECTORY=$1
fi


echo "$DIRECTORY"
icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -g -openmp -I/${GSL_ROOT}/include -I/${MKLROOT}include -I/${MKLROOT}include/fftw -I/${HF5_ROOT}/include -L/${HF5_ROOT}/lib /${HF5_ROOT}/lib/libhdf5_hl_cpp.a /${HF5_ROOT}/lib/libhdf5_cpp.a /${HF5_ROOT}/lib/libhdf5_hl.a /${HF5_ROOT}/lib/libhdf5.a -L/${GSL_ROOT}/lib -L/${MKLROOT}lib/intel64 -Wl,-rpath -Wl,/${HF5_ROOT}/lib -mkl -lgsl -liomp5 -lhdf5 -lhdf5_cpp -lpthread -ldl -lrt -lz -DMKL_Complex16="std::complex<double>"  ${DIRECTORY}/photoion_comp.cpp -o ${DIRECTORY}/PWAPIC_an_orth.exe
#${icpp} -g -openmp  -I/${GSL_ROOT}/include -I/${MKLROOT}include -I/${MKLROOT}include/fftw ${HF5_LINKERS} -L/${GSL_ROOT}/lib -L/${MKLROOT}lib/intel64 -mkl -lgsl -liomp5 -lpthread -ldl -DMKL_Complex16="std::complex<double>"  ${DIRECTORY}/photoion_comp.cpp -o ${DIRECTORY}/photoionization_comp.exe
