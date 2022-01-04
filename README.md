# GraSphInv

## 1 Installation

You need to compile from source. 

### 1.1 Recommended Compiler

The program is written in C++, so a C++ compiler is necessary. You can use one of the following compilers.

- [Intel oneAPI (based on Clang)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.koucg5) (2021, 2022)

- [GNU g++](https://gcc.gnu.org/) (gcc version 9.3.0)

> The versions we have tested are given in brackets

You can modify `Config.cmake` to specify a compiler. The default is oneAPI compiler for which the command is `icpx`.

### 1.2 Other requirements

**3rd party libraries**

- [gsl](https://www.gnu.org/software/gsl/) (necessary)
- [netcdf](https://www.unidata.ucar.edu/software/netcdf/) (optional)

```shell
#The fastest way to install netcdf:

#For ubuntu users
sudo apt-get install libnetcdf-dev
sudo apt-get install libnetcdf-c++4-dev
#For Fedora users
sudo yum install netcdf-devel
sudo yum install netcdf-cxx-devel
```

> For CentOS, old releases of Ubuntu or Fedora, it is recommended to build [netcdf-c](https://www.unidata.ucar.edu/software/netcdf/docs/) and [netcdf-cxx](https://github.com/Unidata/netcdf-cxx4) from the source code. 

**Building tools**

- [GNU make](https://www.gnu.org/software/make/)
- [cmake](https://cmake.org/) 

### 1.3 Build

(1) "cd" the `GraSphInv` directory

(2) Edit `GraSphInv/Config.cmake`

```bash
# Specify the c++ compiler. The default one is icpx (oneAPI).
set(ENV{CXX} icpx)
# Specify whether netcdf library will be used. The value should be TRUE or FALSE
set(USE_NETCDF TRUE)
# Specify whether mkl library will be used. THe value should be TRUE or FALSE
set(USE_MKL TRUE)
```

(2) Type "mkdir build" to create a "build" folder

(3) "cd build" and type cmake ..

(4) Type make

(5) [Optional] add `GraSphInv/build` directory to the environmental variable `PATH` so that you don't need to copy the program files to your working directories.

## 2 Usage



## 3 Output

An inverted model can be  written to three different formats: *.txt, *.vtk, *.nc.

- **txt** Each line represent a cell. But you cannot directly visualized the result in this format.  You need to use python, matlab or GMT to interpolate the values onto regular grids and then plot the results.
- **vtk** The shape of a tesseroid is approximated by a polyhedron in *.vtk file. You can open a vtk file using Paraview. But it is not easy to show axes in spherical coordinates.
- **nc** The netcdf library is used to write results into *.nc files. However, it is difficult to store an irregular mesh into a *.nc file. Here the workaround is to divide the adaptively refined mesh to a regular mesh where the cell size is the smallest one in the adaptively refined mesh



## 4 Examples



