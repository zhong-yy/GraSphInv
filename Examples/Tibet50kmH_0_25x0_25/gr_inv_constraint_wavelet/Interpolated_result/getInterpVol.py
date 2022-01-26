import numpy as np
import interpData as itp

filename='../gr_inv_pet.txt'

# (1) Interpolate the result into a uniform grid with a 15min x 15min x 18km cell size
start_r=5978.137
stop_r=6338.137
num_r=21
start_lon=65
stop_lon=111
num_lon=185
start_lat=20
stop_lat=48
num_lat=113

cols=[6,7,8,10]
skiprows=2
out_file='./model_15m_15m_18km.txt'

LONI, LATI, RI, VI=itp.interp_volume_and_write_out(filename,out_file,
        0.25,0.25,18, 
        start_lon,stop_lon,num_lon,
        start_lat,stop_lat,num_lat,
        start_r,stop_r,num_r,
        cols=cols,skiprows=skiprows,neighbors=50)


# (2) Interpolate the result into a uniform grid with a 30min x 30min x 36km cell size
start_r=5978.137
stop_r=6338.137
num_r=11
start_lon=65
stop_lon=111
num_lon=93
start_lat=20
stop_lat=48
num_lat=57

cols=[6,7,8,10]
skiprows=2
out_file='./model_30m_30m_36km.txt'

LONI, LATI, RI, VI=itp.interp_volume_and_write_out(filename,out_file,
        0.5,0.5,36, 
        start_lon,stop_lon,num_lon,
        start_lat,stop_lat,num_lat,
        start_r,stop_r,num_r,
        cols=cols,skiprows=skiprows,neighbors=50)
