add_library(FilesUtils molpro_parser.cpp) #  hf5_photoion.cpp)

add_library(gsl SHARED IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgsl.so)

add_library(gslcblas SHARED IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgslcblas.so)

find_package (HDF5)
if (HDF5_FOUND)
   target_include_directories(FilesUtils PRIVATE
                                         ${HDF5_INCLUDE_DIRS}
                                         "$ENV{GSL_ROOT}/include"
                                         "${PROJECT_SOURCE_DIR}/Utilities")

   target_link_libraries ( FilesUtils ${HDF5_LIBRARIES}
                                       gsl
                                       gslcblas)
endif (HDF5_FOUND)


#target_link_libraries(FilesUtils PRIVATE gsl)
#target_link_libraries(FilesUtils PRIVATE gslcblas)
#install(TARGETS FilesUtils DESTINATION lib)
#install(FILES files_utils.h DESTINATION include)
