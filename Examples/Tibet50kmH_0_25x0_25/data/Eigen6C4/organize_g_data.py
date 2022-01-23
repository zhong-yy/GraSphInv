import sys
import numpy

Total_data='./GrafLab_output/Tibet_g.txt'
Topo_data='./GrafLab_output/Tibet_topo_g.txt'
Topo_plus_GRS80_data='./GrafLab_output/Tibet_topo_plusGRS80_g.txt'

#read data from files
lats, lons, gx, gy, gz= numpy.loadtxt(Total_data, unpack=True)

lats, lons, gx_topography, gy_topography, gz_topography= numpy.loadtxt(Topo_data, unpack=True)

lats, lons, gx_topography_GRS80, gy_topography_GRS80, gz_topography_GRS80= numpy.loadtxt(Topo_plus_GRS80_data, unpack=True)

gx_normal=gx_topography_GRS80-gx_topography
gy_normal=gy_topography_GRS80-gy_topography
gz_normal=gz_topography_GRS80-gz_topography

f1=open('./Tibet_total_g_r','w')
f2=open('./Normal_g_r','w')
f3=open('./Topo_effect_g_r','w')
f5=open('./Free_air_g_r','w')
f6=open('./Bouguer_g_r','w')

f7=open('./Bouguer_g_theta','w')
f8=open('./Bouguer_g_phi','w')

#free_air=

gx_bouguer=-(gx-gx_topography_GRS80)#x.....-theta
gy_bouguer=-(gy-gy_topography_GRS80)#y.....-phi

gz_freeair=gz-gz_normal
gz_bouguer=gz_freeair-gz_topography#z......r




for i in range(len(lats)):
	f1.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz[i]))
	f2.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_normal[i]))
	f3.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_topography[i]))
	f5.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_freeair[i]))
	f6.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_bouguer[i]))
	
	f7.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gx_bouguer[i]))
	f8.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gy_bouguer[i]))
