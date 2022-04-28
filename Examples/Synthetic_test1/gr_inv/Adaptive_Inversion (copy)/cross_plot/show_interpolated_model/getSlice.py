import interpData as itd
import numpy as np
import matplotlib.pyplot as plt

r0=6238.137
lon0=96.5
lon1=103.5
model_file='../../ada_result.txt'
cols=[6,7,8,10]
skiprows=2
#model_file='./ada_intp_data'
#cols=[0,1,2,4]
#skiprows=0
LAT1, Z1, V1= itd.interp_lon_slice(model_file,lon0,20.0,40.0,81,5978.137,6378.137,41,cols=cols,skiprows=skiprows,neighbors=25,kernel='linear')
LAT_flat1=LAT1.flatten()
Z_flat1=Z1.flatten() 
V_flat1=V1.flatten()
data1=np.hstack((LAT_flat1[:,np.newaxis],Z_flat1[:,np.newaxis],V_flat1[:,np.newaxis]))
print(LAT1.flatten().shape)
print(Z1.flatten().shape)
print(V1.flatten().shape)
np.savetxt('interp_data1.xyz', data1,fmt="%25.15e")

clrmap='rainbow'



LAT2, Z2, V2= itd.interp_lon_slice(model_file,lon1,20.0,40.0,81,5978.137,6378.137,41,cols=cols,skiprows=skiprows,neighbors=25,kernel='linear')
LAT_flat2=LAT2.flatten()
Z_flat2=Z2.flatten() 
V_flat2=V2.flatten()
data2=np.hstack((LAT_flat2[:,np.newaxis],Z_flat2[:,np.newaxis],V_flat2[:,np.newaxis]))
np.savetxt('interp_data2.xyz', data2,fmt="%25.15e")

LON3, LAT3, V3=itd.interp_r_slice(model_file,r0,90.0,110.0,81,20.0,40.0,81,cols=cols,skiprows=skiprows,neighbors=25,kernel='linear')
LON_flat3=LON3.flatten()
LAT_flat3=LAT3.flatten()
V_flat3=V3.flatten()

data3=np.hstack((LON_flat3[:,np.newaxis],LAT_flat3[:,np.newaxis],V_flat3[:,np.newaxis]))

np.savetxt('interp_data3.xyz', data3,fmt="%25.15e")


fig, (ax1,ax2,ax3)=plt.subplots(3,1,figsize=(6,15))
#plt.pcolormesh(LAT1,Z1,V1,cmap=clrmap,shading='gouraud')
sct1=ax1.scatter(LAT1,Z1,c=V1,cmap=clrmap)
fig.colorbar(sct1,ax=ax1,location='right')
sct2=ax2.scatter(LAT2,Z2,c=V2,cmap=clrmap)
fig.colorbar(sct2,ax=ax2,location='right')
sct3=ax3.scatter(LON3,LAT3,c=V3,cmap=clrmap)
fig.colorbar(sct3,ax=ax3,location='right')
ax3.set_aspect('equal', 'box')
plt.show()
