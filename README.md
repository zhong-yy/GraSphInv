# GraSphInv
[![DOI](https://zenodo.org/badge/444264186.svg)](https://zenodo.org/badge/latestdoi/444264186)

## 1 Installation

Compile the program from source. 

### 1.1 Recommended Compiler

The inversion program is written in C++. You can use one of the following compilers.

- [Intel oneAPI (based on Clang)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.koucg5) (2021, 2022)

- [GNU g++](https://gcc.gnu.org/) (gcc version 9.3.0)

> The versions we have tested are given in brackets

The default is oneAPI compiler for which the command is `icpx`. You can modify `Config.cmake` to specify a compiler.

### 1.2 Other requirements

#### 3rd party libraries

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

> Other third-party libraries are headers only, which are included in the `GraSphInv/src/3rd_party_lib` directory and do not need installation.

#### Building tools

- [GNU make](https://www.gnu.org/software/make/)
- [cmake](https://cmake.org/) 

### 1.3 Build

> NOTE: If you don't use the oneAPI c++ compiler, or you don't want to install netcdf, you need to modify  `GraSphInv/Config.cmake`  before compiling.

(1) `cd` the `GraSphInv` directory

(2) `mkdir build`

(3) `cd build`

(4) `cmake ..`

(5) `make`

(6) `make install`

(7) [*Optional*] add `GraSphInv/bin` directory to the environmental variable `PATH` in `~/.bashrc`, so that you don't need to copy the  executable files to your working directories every time.

## 2 Usage

### Inversion

After compiling the source code, you will see a program called GraSphInv. It is the program for inversion

```
GraSphInv [configuration file name]
```

The main configuration file is like this:

```
config_data
config_model
config_inversion
```

where config_data is a file that describes input data,  config_inversion is file that specifies inversion parameters, config_model is a file that specify the inversion region and initial mesh discretization.

#### Some examples

- See `Examples/Synthetic_test1/gr_inv/Adaptive_Inversion` as an example for inversion using the an adaptively refined invesion mesh.

- See `Examples/Synthetic_test1/gr_inv/Adaptive_Inversion/Adaptive_Inversion_crg` and `Examples/Synthetic_test1/gr_inv/Adaptive_Inversion/Adaptive_Inversion_pet` as an example using constraints

- See `Examples/Synthetic_test1/gr_inv/Non-adaptive_Inversion`as an example for inversion using an fixed uniform mesh.
- See `Examples/Synthetic_test1/gr_inv_wavelet` or `Examples/Synthetic_test1/gr_inv_wavelet/Non-adaptive_Inversion` as an example for inversion using wavelet compression

### Crustal correction

RemoveCrustalEffect is a program for crustal correction. See `Examples/Tibet50kmH_0_25x0_25/data/Crust_Correction` as an example.

## 4 Output

An inverted model can be  written to three different formats: *.txt, *.vtk, *.nc.

- **txt** Each line represent a cell. But you cannot directly visualized the result in this format.  You need to use python, matlab or GMT to interpolate the values onto regular grids and then plot the results.
- **vtk** The shape of a tesseroid is approximated by a polyhedron in *.vtk file. You can open a vtk file using Paraview. But it is not easy to show axes in spherical coordinates.
- **nc** The netcdf library is used to write results into *.nc files. However, it is difficult to store an irregular mesh into a *.nc file. Here the workaround is to divide the adaptively refined mesh to a regular mesh where the cell size is the smallest one in the adaptively refined mesh

## 5 References

Yiyuan Zhong, Zhengyong Ren, Jingtian Tang, Yufeng Lin, Bo Chen, Yangfan Deng. "Constrained gravity inversion with adaptive inversion grid refinement in spherical coordinates and its application for  mantle structure beneath Tibetan plateau", under review, 2022



