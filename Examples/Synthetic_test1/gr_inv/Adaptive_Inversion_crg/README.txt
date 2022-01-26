To run the inversion, type in the terminal:

GraSphInv config


File tree structrue in this directory:
.
├── ada_result_crg.nc  :inverted model in netcdf file format.
├── ada_result_crg.txt :inverted model in text format
├── ada_result_crg.vtk :inverted model in vtk format
├── config             :configuration file
├── config_data        :configuration file for observed data
├── config_inversion   :configuration file for inversion parameters
├── config_model       :configuration file for the initial inversion mesh
├── crg_model          :the velocity model given in text format. Each line represents: latitude longitude depth (km) value
├── cross_plot
├── dobs_g_r           :observed data
├── dpredicted_g_r     :predicted data
├── info
├── initial_mesh.vtk
├── Iteration_misfit_GN
├── lambda_misfit_GN
├── log
├── plot_misfit.py
├── plots              :plot the inversion result
└── README.txt
