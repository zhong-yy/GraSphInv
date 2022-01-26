3 directories
.
├── gr_inv                   : Inversion using full sensitivity matrices
├── gr_inv_wavelet           : Inversion using wavelet compression
└── True_model               : plot the synthetic model


8 files
.
├── convert_ref.py           : convert the S-wave model to density model
├── convert_test.py
├── crg_model                : S-wave velocity model
├── dobs_exact_g_r           : Synthetic data without noice
├── dobs_g_r                 : Synthetic data with random error of which the standard variance is 2%*d+0.01*max(d) (d means data values)
├── README.txt
├── ref_model                : density model converted from the S-wave model using the python script 'convert_ref.py'
└── sites                    : observation points
