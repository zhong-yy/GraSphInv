import sys
import numpy

Total_data='Tibet_g_ggt.txt'
Topo_data='Tibet_topo_g_ggt.txt'
Topo_plus_GRS80_data='Tibet_topo_plusGRS80_g_ggt.txt'

#read data from files
lats, lons, gx, gy, gz, Vxx, Vyy, Vzz, Vxy, Vxz, Vyz = numpy.loadtxt(Total_data, unpack=True)

lats, lons, gx_topography, gy_topography, gz_topography, Vxx_topography, Vyy_topography, Vzz_topography, Vxy_topography, Vxz_topography, Vyz_topography = numpy.loadtxt(Topo_data, unpack=True)

lats, lons, gx_topography_GRS80, gy_topography_GRS80, gz_topography_GRS80, Vxx_topography_GRS80, Vyy_topography_GRS80, Vzz_topography_GRS80, Vxy_topography_GRS80, Vxz_topography_GRS80, Vyz_topography_GRS80 = numpy.loadtxt(Topo_plus_GRS80_data, unpack=True)

gx_normal=gx_topography_GRS80-gx_topography
gy_normal=gy_topography_GRS80-gy_topography
gz_normal=gz_topography_GRS80-gz_topography
Vzz_normal=Vzz_topography_GRS80-Vzz_topography
Vxx_normal=Vxx_topography_GRS80-Vxx_topography
Vyy_normal=Vyy_topography_GRS80-Vyy_topography
Vxy_normal=Vxy_topography_GRS80-Vxy_topography
Vxz_normal=Vxz_topography_GRS80-Vxz_topography
Vyz_normal=Vyz_topography_GRS80-Vyz_topography

f1=open('./Tibet_total_g_r','w')
f2=open('./Normal_g_r','w')
f3=open('./Topo_effect_g_r','w')
f5=open('./Free_air_g_r','w')
f6=open('./Bouguer_g_r','w')

f8=open('./Tibet_total_T_rr','w')
f9=open('./Normal_T_rr','w')
f10=open('./Topo_effect_T_rr','w')
f12=open('./Free_air_T_rr','w')
f13=open('./Bouguer_T_rr','w')
#f15=open('./Bouguer_gr_Trr','w')
f16=open('./Bouguer_T_theta_theta','w')
f17=open('./Bouguer_T_phi_phi','w')
f18=open('./Bouguer_T_theta_phi','w')
f19=open('./Bouguer_T_theta_r','w')
f20=open('./Bouguer_T_phi_r','w')

f21=open('./Bouguer_g_theta','w')
f22=open('./Bouguer_g_phi','w')

#free_air=

gx_bouguer=-(gx-gx_topography_GRS80)#x.....-theta
gy_bouguer=-(gy-gy_topography_GRS80)#y.....-phi

gz_freeair=gz-gz_normal
gz_bouguer=gz_freeair-gz_topography#z......r

Vzz_freeair=Vzz-Vzz_normal
Vzz_bouguer=Vzz_freeair-Vzz_topography#z......r
Vxx_bouguer=Vxx-Vxx_topography_GRS80
Vyy_bouguer=Vyy-Vyy_topography_GRS80
Vxy_bouguer=Vxy-Vxy_topography_GRS80
Vxz_bouguer=-(Vxz-Vxz_topography_GRS80)
Vyz_bouguer=-(Vyz-Vyz_topography_GRS80)




for i in range(len(lats)):
	f1.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz[i]))
	f2.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_normal[i]))
	f3.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_topography[i]))
	f5.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_freeair[i]))
	f6.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gz_bouguer[i]))

	f8.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vzz[i]))
	f9.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vzz_normal[i]))
	f10.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vzz_topography[i]))
	f12.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vzz_freeair[i]))
	f13.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vzz_bouguer[i]))
#	f15.write(" %.12e  %.12e %.12e %.12e \r\n" % (lats[i], lons[i],gz_bouguer[i], Vzz_bouguer[i]))
	
	f16.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vxx_bouguer[i]))
	f17.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vyy_bouguer[i]))
	f18.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vxy_bouguer[i]))
	f19.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vxz_bouguer[i]))
	f20.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], Vyz_bouguer[i]))
	
	f21.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gx_bouguer[i]))
	f22.write(" %.12e  %.12e %.12e \r\n" % (lats[i], lons[i], gy_bouguer[i]))
