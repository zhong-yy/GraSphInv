import sys
import interpData as itd
import numpy as np

filename=sys.argv[1]
interp_coords_file=sys.argv[2]
output_filename=sys.argv[3]

cols=[6,7,8,10]
itd.interp_profiles(filename,interp_coords_file,output_filename,cols=cols,skiprows=2,neighbors=100,kernel='linear',ref_surface=6378.137)


