import sys
import interpData as itd
import numpy as np
#def interp_r_slice(filename,r0,start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
filename=sys.argv[1]
r0=float(sys.argv[2])
start_lon=float(sys.argv[3])#65
stop_lon=float(sys.argv[4])#111
num_lon=int(sys.argv[5]) #165
start_lat=float(sys.argv[6]) #20
stop_lat=float(sys.argv[7]) #48
num_lat=int(sys.argv[8]) #123
cols=[6,7,8,10]
LON,LAT,V=itd.interp_r_slice(filename,r0,start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,cols=cols,skiprows=2,neighbors=100,kernel='linear')
#thin_plate_spline

LON_flat=LON.flatten()
LAT_flat=LAT.flatten()
V_flat=V.flatten()

data=np.hstack((LON_flat[:,np.newaxis],LAT_flat[:,np.newaxis],V_flat[:,np.newaxis]))
output_filename=sys.argv[9]
np.savetxt(output_filename, data, fmt="%25.15e")

