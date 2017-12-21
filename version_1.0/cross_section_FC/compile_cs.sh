#!/bin/bash
icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
${icpp} -qopenmp  -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${MKLROOT}lib/intel64 -mkl -liomp5 -lpthread -ldl  /data1/home/stephan/photoionization_coupling_comp/version_1.0/cross_section_FC/cross_section_FC.cpp -o /data1/home/stephan/photoionization_coupling_comp/version_1.0/cross_section_FC/total_PICE_comp.exe
