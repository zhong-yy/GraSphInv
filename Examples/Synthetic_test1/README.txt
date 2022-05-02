5 directories
.
├── ggt_inv                  : Inversion of gravity gradient tensor
├── gr_inv                   : Inversion of gr data using wavelet compression
├── gr_inv_wavelet           : Inversion of gr data using full sensitivity matrices
├── Trr_inv                  : Inversion of Trr data
└── True_model               : plot the synthetic model








files
.
├── convert_ref.py           : convert the S-wave model to density model
├── convert_test.py
├── crg_model                : S-wave velocity model
├── dobs_exact_g_r           : Synthetic data without noice
├── dobs_g_r                 : Synthetic data with random error of which the standard variance is 2%*d+0.01*max(d) (d means data values)
├── README.txt
├── ref_model                : density model converted from the S-wave model using the python script 'convert_ref.py'
└── sites                    : observation points
