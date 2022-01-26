.
├── dobs_exact_g_r   :this file is the observed data without noise. It is only used to show the distribution of the added noise
├── get_mean_std.py  :a python script calculating means and stds
├── make_frame.py    :generate a closed "box" as an input of 'gmt plot' to show the outline of the true model
├── max_min.py
├── plot_comparison.sh    :plot a comparison beween the results from the adaptive mesh and the fixed uniform mesh
├── plot_field2.sh        :plot the observed data and the predicted data
├── plot_noise.sh         :show the noise distribution in the observed data
├── plot_rho_same_colorbar.sh   :plot the inversion result in three cross-sections using the same colorbar
├── plot_rho.sh                 :plot the inversion result in three cross-sections. Each cross-section uses an independent colorbar
├── README.txt
├── swap_2cols_of_g_data.py     :make format of the observed/predicted gravity data suitable for GMT (e.g. column 1: longitude, column 2: latitude)
├── synthetic_test_gr_misfit.jpg
├── synthetic_test_gr_model_no_cst.eps
└── synthetic_test_gr_model_no_cst.jpg




1 Open a terminal in this directory, and type

bash plot_field2.sh

to plot the data.

2 Type 

bash plot_rho_same_colorbar.sh

to plot the inversion result.

3 Type 

bash plot_comparison.sh

to compare the inverted models from the adaptive mesh and the fixed mesh in the same figure.
