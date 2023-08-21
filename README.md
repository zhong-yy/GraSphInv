# GraSphInv

[![DOI](https://zenodo.org/badge/444264186.svg)](https://zenodo.org/badge/latestdoi/444264186)

**GraSphInv** is a program for gravity inversion using adaptive inversion mesh refinement in the spherical coordinate system.  A-priori constraints can be used in the inversion. Wavelet compression can be used to reduce memory consumption. 

The input data can be any combination of components of vector gravity field or gravity gradient tensor. The output is a 3D distribution of density contrasts. 

## 1 Installation

### 1.1 Prerequisites

#### 1.1.1 Compiler

The inversion program is written in C++. One of the following compilers can be used to compile the program:

- [Intel oneAPI (based on Clang)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.koucg5) (2021, 2022 **recommended**)
- [GNU g++](https://gcc.gnu.org/) (gcc version 9.3.0, not fully tested)

> The versions we have tested are given in brackets

The default is oneAPI compiler for which the command is `icpx`. You can modify the file `GraSphInv/Config.cmake` to specify a compiler. For example, change `set(ENV{CXX} icpx)` to `set(ENV{CXX} g++)` to use a `g++` compiler. Since MKL is only supported by intel compilers, if the `g++` compiler is used, you may need to change `set(USE_MKL TRUE)` to `set(USE_MKL FALSE)` in file `GraSphInv/Config.cmake`.

#### 1.1.2 Third-party libraries

(1) Eigen

We use the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for basic linear algebra operations. Eigen is a template library, so you don't need to worry about the installation. A copy of Eigen can be found in `GraSphInv/src/3rd_party_lib/eigen`.

(2) GNU Scientific Library (GSL)

**[gsl](https://www.gnu.org/software/gsl/)** is only used to calculate wavelet transforms. Follow the instructions in the gsl source code package to install it (some usual installation steps, `configure`, `make`, `make install`).

(3) NetCDF (optional)

Netcdf library is used to write the inverison model to netcdf4 format (`.nc` files). But if you don't need netcdf4 outputs, it is not necessary to compile this inversion program with `netcdf` library. 

> A downside of netcdf format is that before a final model is written into a netcdf file, all cells in the irregular model mesh need to be subdivided into the size of the smallest cell in the mesh. However, the irregular mesh is kept as it is in the `.vtk` or `.txt` outputs.

```shell
#The fastest way to install netcdf:
#For ubuntu users
sudo apt-get install libnetcdf-dev
sudo apt-get install libnetcdf-c++4-dev
#For Fedora users
sudo yum install netcdf-devel
sudo yum install netcdf-cxx-devel
#For CentOS, old releases of Ubuntu or Fedora, 
# it is recommended to build netcdf-c netcdf-cxx from the source code.
# https://www.unidata.ucar.edu/software/netcdf/docs/
# https://github.com/Unidata/netcdf-cxx4
```

To disable netcdf library before building the program,  open the file `GraSphInv/Config.cmake`, change 

```
set(USE_NETCDF TRUE)
```

to

```
set(USE_NETCDF FALSE)
```

### 1.2 Build

The building tools [GNU make](https://www.gnu.org/software/make/) and [cmake](https://cmake.org/) should be installed before building the program, which are available on most linux platforms. 

Steps to build the program:

(1) `cd` the `GraSphInv` directory

(2) `mkdir build`

(3) `cd build`

(4) `cmake ..`

(5) `make`

(6) `make install`

After `make install`, all executable programs can be found in the `GraSphInv/bin` directory. Optionally, one may want to add  `GraSphInv/bin` directory to the environmental variable `PATH` in `~/.bashrc`, so that  it is not necessary to copy the  executable files to working directories every time. For example, open the `~/.bashrc` file

```bash
$ vim ~/.bashrc
```

and add the following lines:

```bash
export PATH=/home/yyzhong/GraSphInv/bin${PATH:+:${PATH}}                                                                                                                                                           
export PYTHONPATH=[path_to_GraSphInv]/GraSphInv/bin${PYTHONPATH:+:${PYTHONPATH}}
```

where `[path_to_GraSphInv]` should be changed the path to `GraSphInv` in your computer. You can open a terminal under the GraSphInv folder and type `pwd` to check the path.

## 2 Usage

### 2.1 Inversion

After compiling the source code, you will see a program called GraSphInv. The command is the program name followed by the filename of a configuration file.

```
GraSphInv [configuration_file_name]
```

The configuration file contains only three lines: the first line is the filename of a configuration file for data, the second line is the filename of the configuration file for the inversion region, the third line is the filename of the configuration file for inversion parameters. Then a main configuration file is like this:

```
config_data
config_model
config_inversion
```

where **config_data** is the file that describes input data,  **config_inversion** is the file for inversion parameters,  and **config_model** is the file that specify the inversion region and initial mesh discretization.

### 2.2 Crustal correction

RemoveCrustalEffect is a program for crustal correction. See `Examples/Tibet50kmH_0_25x0_25/data/Crust_Correction` as an example.

## 3. Some examples

### 3.1 Inversion of gr data

(1) Preparation of synthetic data

There two executable programs for generating synthetic observed data. `Synthetic_data` is used to generate gravity anomalies, and `Synthetic_data_ggt` is used to generate gravity gradient tensor data.

```bash
$ cd Examples/Synthetic_test1
$ Synthetic_data #generate gravity data
$ Synthetic_data_ggt #generate gravity gradient data
```

(2) Preparation of synthetic constaint models.

After running `Synthetic_data`, you can see a data file **dobs_g_r**, and a model file of a velocity model used as constraint: **crg_model**. Use the python script to convert the velocity model to a reference density model: 

```bash
python convert_ref.py
```

Now, you have a data file **dobs_g_r**, a reference density model file **ref_model**, a velocity constraint model file **crg_model**.

(3) Inversion using gr component, without a-priori constraints.

```bash
# copy the data file to a separate directory
$ cp ./dobs_g_r ./gr_inv/Adaptive_Inversion
# change to that directory
$ cd ./gr_inv/Adaptive_Inversion
```

You will see four configuration files in the **Adaptive_Inversion** folder:

```
├── config              : configuration file
├── config_data         : data configuration file
├── config_inversion    : configuration for inversion parameters
└── config_model        : configuration for the inversion region
```

There are comments in these file, you can open them to see the formats and usage of them.

Next, run the inversion:

```bash
$ GraSphInv config
```

> If GraSphInv/bin has not been added to the environment variable PATH, then you may need to copy it to the present working directory.

The resulting model will be stored in ada_result.txt, ada_result.vtk and ada_result.nc files. 

Then plot the inversion result

```bash
$ cd plots
# show the inverted model
$ bash plot_rho_same_colorbar.sh
# or
$ bash plot_rho.sh
# show the comparison between the observed data and the predicted data
$ bash plot_field.sh
```

In file `gr_inversion/Adaptive_Inversion/config_inversion`, line 50 specify whether the mesh will be refined:

```bash
#maximum times of refinement. If it's 0, the mesh will not be refined.
8
```

You can compare it with `Non-adaptive_Inversion/config_inversion` and see the difference.

> See `Examples/Synthetic_test1/gr_inv/Non-adaptive_Inversion`as an example for inversion using an fixed uniform mesh.

### 3.2 Inversion with a-priori constraints

(1) Using a reference and initial model (which may be from direct paramter relationshisp between density and velocity)

```bash
# go to the Synthetic_test1 folder
$ cd Examples/Synthetic_test1
# copy the data file
$ cp ./dobs_g_r ./gr_inv/Adaptive_Inversion_pet
# copy the reference model file
$ cp ./ref_model ./gr_inv/Adaptive_Inversion_pet
# go to the directory Adaptive_Inversion_pet
$ cd ./gr_inv/Adaptive_Inversion_pet

# have a look at the configuration files
$ vim config
$ vim config_data
$ vim config_model
$ vim config_inversion

# run the inversion
$ GraSphInv config
```

Line 76 in file `config_inversion`

```bash
ref_model xyz 40 40 20 0 
```

specifies the reference model file **ref_model**.  `xyz` means the coordinates are ordered as `latitude, longitude, depth`. If the order is  "radius, latitude, longitde", then it should be `zxy`. In addition, `40 40 20` are the number of grid nodes along latitude, longitude and depth, whose order is also subject to the order of `xyz`. The last number represents the coordinate that changes fastest in the file, with 0 indicating latitude changes fastest, 1 indicating longitude changes fastest.

For example, the following file should be described by `xyz 3 2 1 0`

```bash
# lat lon dep value
20   100   50
21   100   50
22   100   50
20   101   50
21   101   50
22   101   50
```

The following file should be described by `yxz 2 3 1 1`

```bash
# lon lat dep
100 20 50
101 20 50
100 21 50
101 21 50
100 22 50
101 22 50
```

(2) Cross gradient constraint

```bash
# go to the Synthetic_test1 folder
$ cd Examples/Synthetic_test1
# copy the data file
$ cp ./dobs_g_r ./gr_inv/Adaptive_Inversion_crg
# copy the velocity model file
$ cp ./crg_model ./gr_inv/Adaptive_Inversion_crg
# go to the directory Adaptive_Inversion_crg 
$ cd ./gr_inv/Adaptive_Inversion_crg

# have a look at the configuration files
$ vim config
$ vim config_data
$ vim config_model
$ vim config_inversion

# run the inversion
$ GraSphInv config
```

### 3.3 Inversion with wavelet compression

See `Examples/Synthetic_test1/gr_inv_wavelet` or `Examples/Synthetic_test1/gr_inv_wavelet/Non-adaptive_Inversion` as an example for inversion using wavelet compression.

Lines 13 and 16 in `gr_inv_wavelet/Adaptive_Inversion` specify whether to use wavelet compression and a relative theshold which affects the compression ratio. A larger threshold leads to more compression, and also larger approximation error.

```
# 0 or 1: full sensitivity matrix
# 2: wavelet compression
2
0.005
```

### 3.4 Inversion of gravity gradient tensor

Go to the Synthetic_test1 directory and run

```bash
$ Synthetic_data_ggt
```

to get synthetic gravity gradient tensor data.

Copy the GGT data to `ggt_inv` directory.

```bash
$ cp ./dobs_T_rr ./dobs_T_rphi ./dobs_T_rtheta ./dobs_T_thetatheta ./dobs_T_thetaphi ./dobs_T_phiphi ./ggt_inv
```

```bash
$ cd ./ggt_inv
```

Combine 5 GGT components into one file, 

```bash
python combine5files.py
```

Specify the number of used components, and the order of different components in the data configuration file `config_data`:

```
...
#How many components will be use?
5

#gravity component markers
#0:V   1:gr   2:g_theta   3:g_phi   4:T_rr   5:T_rtheta   6:T_rphi   
#7:T_thetatheta   8:T_thetaphi   9:T_phiphi
4 5 6 7 8
...
```

Run the inversion

```bash
GraSPhInv config
```

Then use the shell scripts to in the `/plot/` subdirectory to display the inversion results.

## 4 Output

### 4.1 Description

An inverted model can be  written to three different formats: *.txt, *.vtk, *.nc.

- **txt** Each line represent a cell. But you cannot directly visualized the result in this format.  You need to use python, matlab or GMT to interpolate the values onto regular grids and then plot the results.
- **vtk** The shape of a tesseroid is approximated by a polyhedron in *.vtk file. You can open a vtk file using Paraview, but it is triky to view spherical coordinates in paraview.
- **nc** The netcdf library is used to write results into *.nc files. However, it is difficult to store an irregular mesh into a *.nc file. Here the workaround is to divide the adaptively refined mesh into a regular mesh where each cell has the same size as the smallest cell in the adaptively refined mesh.

### 4.2 Visualization

GMT, python matplotlib, Paraview can be used to visualize the inversion results.

We also provide a python script `interpData.py` to interpolate inverted values at 3D scattered points to regular meshes on 2D cross-sections. To use `interpData.py`, scipy >=1.7.2 is required.

## 5 References

Yiyuan Zhong, Zhengyong Ren, Jingtian Tang, Yufeng Lin, Bo Chen & Yangfan Deng, Yingde Jiang (2022). Constrained gravity inversion with adaptive  inversion grid refinement in spherical coordinates and its application  to mantle structure beneath Tibetan Plateau. *Journal of Geophysical Research: Solid Earth*, 127, e2021JB022916. https://doi.org/10.1029/2021JB022916

```tex
@article{zhong2022,
author = {Zhong, Yiyuan and Ren, Zhengyong and Tang, Jingtian and Lin, Yufeng and Chen, Bo and Deng, Yangfan and Jiang, Yingde},
title = {Constrained Gravity Inversion With Adaptive Inversion Grid Refinement in Spherical Coordinates and Its Application to Mantle Structure Beneath Tibetan Plateau},
journal = {Journal of Geophysical Research: Solid Earth},
volume = {127},
number = {5},
pages = {e2021JB022916},
doi = {https://doi.org/10.1029/2021JB022916},
year = {2022}
}
```

## Troubleshooting
1. Error while loading shared libraries XXX.so.XX cannot open shared object file: No such file or directory

Try adding`export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH` to ~/.bashrc

2. Unexpected computation errors in the `Tibet50kmH_0_25x0_25/data/Crust_Correction` example

If you are using g++ 11.3 or 9.5 and encounter this issue, please try intel compiler first.  See issue https://github.com/zhong-yy/GraSphInv/issues/1. I still have no idea of the exact cause of this problem. 

Intel oneAPI 2021,2022 have been tested. Although I don't have time to make a comprehensive test of Intel oneAPI 2023 on all examples, I have tested it on the crust correction example and I suppose it also works well for other examples. However, efficacy of future version of Intel compiler is not guaranteed.

