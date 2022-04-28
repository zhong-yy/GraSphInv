import numpy as np
import interpData as itp


filename='../ada_result.txt'
start_r=5988.137
stop_r=6368.137
num_r=20

start_lon=90.25
stop_lon=109.75
num_lon=40

start_lat=20.25
stop_lat=39.75
num_lat=40

cols=[6,7,8,10]
skiprows=2
out_file='./intp_vol_data_ada.txt'
LONI, LATI, RI, VI=itp.interp_volume_and_write_out(filename,out_file,
        0.5,0.5,20, 
        start_lon,stop_lon,num_lon,
        start_lat,stop_lat,num_lat,
        start_r,stop_r,num_r,
        cols=cols,skiprows=skiprows,neighbors=50)

## Result of the inversion using a uniform mesh
filename2='../../Non-adaptive_Inversion/non_ada_result.txt'

out_file2='./intp_vol_data_non_ada.txt'

LONI2, LATI2, RI2, VI2=itp.interp_volume_and_write_out(filename2,out_file2,0.5,0.5,20, start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,start_r,stop_r,num_r,cols=cols,skiprows=skiprows,neighbors=50)
