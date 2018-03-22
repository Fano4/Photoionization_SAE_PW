#!/bin/bash
#icpp=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
icpp=icpc
#${icpp} -qopenmp  -I/${MKLROOT}include -I/${MKLROOT}include/fftw  -L/${MKLROOT}lib/intel64 -mkl -liomp5 -lpthread -ldl  /CECI/home/ulg/cpt/svdwild/photoionization_v1.0/cros_section_integration_cube.cpp -o /CECI/home/ulg/cpt/svdwild/photoionization_v1.0/cs_computation.exe
${icpp}  /CECI/home/ulg/cpt/svdwild/photoionization_v1.0/cros_section_integration_cube.cpp -o /CECI/home/ulg/cpt/svdwild/photoionization_v1.0/cs_computation.exe
