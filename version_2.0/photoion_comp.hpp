//standard libraries
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include "mkl.h"
#include "omp.h"

//internal program files
#include "objects.hpp"
#include "files_reader.hpp"
#include "objects.cpp"
#include "files_reader.cpp"

//functions_declaration
int photoion_comp(int argc,char* argv[]);
