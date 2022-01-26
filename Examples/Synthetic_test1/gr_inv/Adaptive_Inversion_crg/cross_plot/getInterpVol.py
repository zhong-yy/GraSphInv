import numpy as np
import interpData as itp


filename='../ada_result_crg.txt'
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
#LONI, LATI, RI, VI=itp.interp_volume(filename, start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,start_r,stop_r,num_r,cols=cols,skiprows=skiprows,neighbors=50)

#lon_flat=LONI.flatten()
#lat_flat=LATI.flatten()
#r_flat=RI.flatten()
#V_flat=VI.flatten()
#
#data=np.hstack((lon_flat[:,np.newaxis],lat_flat[:,np.newaxis],r_flat[:,np.newaxis],V_flat[:,np.newaxis]))
#
#np.savetxt(out_file,data,fmt="%25.15e")


## Result of the inversion using a uniform mesh
filename2='../../Non-adaptive_Inversion_crg/non_ada_result_crg.txt'

out_file2='./intp_vol_data_non_ada.txt'

LONI2, LATI2, RI2, VI2=itp.interp_volume_and_write_out(filename2,out_file2,0.5,0.5,20, start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,start_r,stop_r,num_r,cols=cols,skiprows=skiprows,neighbors=50)

#lon_flat2=LONI2.flatten()
#lat_flat2=LATI2.flatten()
#r_flat2=RI2.flatten()
#V_flat2=VI2.flatten()
#
#data2=np.hstack((lon_flat2[:,np.newaxis],lat_flat2[:,np.newaxis],r_flat2[:,np.newaxis],V_flat2[:,np.newaxis]))
#
#np.savetxt(out_file2,data2,fmt="%25.15e")
