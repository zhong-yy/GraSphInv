.
├── ada_result.nc    :inverted model in netcdf file format.
├── ada_result.txt   :inverted model in text format
├── ada_result.vtk   :inverted model in vtk format
├── config           :configuration file
├── config_data      :configuration file for observed data
├── config_inversion :configuration file for inversion parameters
├── config_model     :configuration file for the initial inversion mesh
├── cross_plot       :make a crossplot of the inverted values from the inversion mesh and the fixed mesh
├── dobs_g_r         :observed data
├── dpredicted_g_r   :predicted data from the inverted model
├── info             :a record of the performance of the inversion
├── initial_mesh.vtk     :initial mesh written in vtk format
├── Iteration_misfit_GN  
├── lambda_misfit.eps
├── lambda_misfit_GN
├── log
├── max_min.py
├── plot_misfit.py      :plot the data misfit vs. lambda values
├── plots               :make figures of the inversion result
├── README.txt
├── refinement_process  :visulization of the refinement process
├── result_at_0.nc      :result of the first iteration
├── result_at_0.vtk     :result of the first iteration
├── result_at_1.nc      :result of the second iteration
├── result_at_1.vtk     :result of the second iteration 
├── result_at_2.nc      :result of the third iteration 
└── result_at_2.vtk     :result of the third iteration 


To run the inversion, type in the terminal:

GraSphInv config
