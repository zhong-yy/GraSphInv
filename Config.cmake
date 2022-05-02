# Specify the c++ compiler. The default one is icpx (oneAPI).
#set(ENV{CXX} g++)
set(ENV{CXX} icpx)

# Specify whether netcdf library will be used. The value should be TRUE or FALSE
set(USE_NETCDF TRUE)

# Specify whether mkl library will be used. THe value should be TRUE or FALSE
set(USE_MKL TRUE)

#set (CMAKE_INSTALL_PREFIX "prefix_path")
