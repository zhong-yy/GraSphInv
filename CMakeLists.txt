

# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --versionh
cmake_minimum_required(VERSION 3.10)

include("Config.cmake")

# set the project name
project(GraSphInv
    VERSION 1.0.0)

set (CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})


set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
find_package(OpenMP REQUIRED)

message("CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}")
message("Compiler command: $ENV{CXX}")
message("Openmp flag ${OpenMP_CXX_FLAGS}")
message("Installation location: ${CMAKE_INSTALL_PREFIX}/bin")
if (${USE_NETCDF})
    message("USE_NETCDF=${USE_NETCDF}; The netcdf library will be included.")
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
#message("${CMAKE_CXX_COMPILER}")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()


aux_source_directory(src/Geometry MESH_SRCS)
aux_source_directory(src/Inversion INVERSION_SRCS)
add_library(Inversion_objs OBJECT
    ${MESH_SRCS}
    ${INVERSION_SRCS}
    )

target_include_directories(Inversion_objs
    PUBLIC
        src/Geometry
        src/Inversion
        src/3rd_party_lib/boost_1_72_0
        src/3rd_party_lib/eigen
        src/3rd_party_lib/Linterp
        src/3rd_party_lib/Timer
        #if netcdf library is used
        $<$<BOOL:$ENV{USE_NETCDF}>:
            src/3rd_party_lib/netcdf/netcdf-cxx4-4.3.1
            src/3rd_party_lib/netcdf/netcdf-c-4.7.4
        >
        )
target_compile_options(Inversion_objs
    PUBLIC
    -Wfatal-errors 
    ${OpenMP_CXX_FLAGS}
    $<$<BOOL:${USE_MKL}>:-qmkl=parallel>
    )
target_compile_definitions(Inversion_objs
    PUBLIC
    $<$<CXX_COMPILER_ID:Intel>:USE_MKL>
    $<$<BOOL:${USE_NETCDF}>:USE_NETCDF>
    )
target_link_libraries(Inversion_objs
        PUBLIC
        OpenMP::OpenMP_CXX
        $<$<BOOL:${USE_NETCDF}>:netcdf netcdf_c++4>
        gsl
        gslcblas
        )
target_link_options(Inversion_objs
    PUBLIC
    $<$<BOOL:${USE_MKL}>:-qmkl=parallel -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl>
    )
##################################executables##########################################
add_executable(GraSphInv
    src/programs/GraSphInv.cpp)
target_link_libraries(GraSphInv
    PRIVATE
    Inversion_objs
    )

add_executable(Synthetic_data
    src/programs/Synthetic_data.cpp)
target_link_libraries(Synthetic_data
    PRIVATE
    Inversion_objs
    )
add_executable(Synthetic_data_ggt
    src/programs/Synthetic_data_ggt.cpp)
target_link_libraries(Synthetic_data_ggt
    PRIVATE
    Inversion_objs
    )

add_executable(RemoveCrustalEffect
    src/programs/RemoveCrustalEffect.cpp)
target_link_libraries(RemoveCrustalEffect
    PRIVATE
    Inversion_objs
    )

add_executable(ConvertTxt2Vtk
    src/programs/ConvertTxt2Vtk.cpp)
target_link_libraries(ConvertTxt2Vtk
    PRIVATE
    Inversion_objs
    )

add_executable(Extract_slice
    src/programs/Extract_netcdf_slice.cpp)
target_link_libraries(Extract_slice
    PRIVATE
    Inversion_objs
    )

add_executable(Test_ordering
    src/programs/Test_ordering.cpp)
target_link_libraries(Test_ordering
    PRIVATE
    Inversion_objs
    )

add_executable(test_wavelet_forward
    src/programs/test_wavelet_forward.cpp)
target_link_libraries(test_wavelet_forward
    PRIVATE
    Inversion_objs
    )
install(FILES py_src/interpData.py
    TYPE BIN)
install(TARGETS GraSphInv RemoveCrustalEffect Synthetic_data Synthetic_data_ggt ConvertTxt2Vtk Extract_slice Test_ordering test_wavelet_forward
    )
#if(UNIX)
#message("Add the installation path to your PATH:\n1. Edit ~/.bashrc file, and add the following line\nexport PATH=${CMAKE_INSTALL_PREFIX}/bin$\{PATH:+:$\{PATH}}")
#endif()
