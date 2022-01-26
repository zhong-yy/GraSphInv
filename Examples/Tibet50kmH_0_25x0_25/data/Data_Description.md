1. **Gravity data**

Gravity data extracted from EIGEN-6C4 are in the folder Eigen6C4/GrafLab_output

2. **Data processing**

(1) Bouguer anomaly

Eigen6C4 contains the Bouguer gravity anomaly, which is processed by the python script Eigen6C4/organize_g_data.py

(2) Crustal correction

The folder **Crust_Correction** contains corrected anomalies.  cal.sh is a script that runs the correction program. In the file cal.sh, `RemoveCrustalEffect` is the name of the binary executable 

- mantle_g_r is the gravity anomaly after removing all effects of sedimentary layers, crystalline crust layers, and Moho. 

- sediments_g_r, crystalline_crust_g_r, moho_g_r  are gravity effects of sediments, crystalline crust and Moho undulation. 

- model_visualization_vtk.tar.xz contains the models for visualization. But Earth's radius is far more large than layer thickness, so that the layer is almost invisible after slicing the model. You can zoom in to see these layers.

(3) Compare processed anomalies from EIGEN-6C4 and EGM2008

Enter the folder **CompareEGM2008**, run the script EIGEN_EGM2.sh:

```
bash EIGEN_EGM2.sh
```



