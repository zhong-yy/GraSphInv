#plot interpolated data
import numpy as np
import sys
#import matplotlib.pyplot as plt
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
#import matplotlib.colors as mcolors
#import matplotlib.gridspec as gridspec
#
##own module
#import interpData
#import similarity

#comparison_method=similarity.ssim
#comparison_method=similarity.normalized_cross_correlation
#win_size=(10,15)
def cross_gradient(X,Z,m,s):
#    nz,nx=X.shape
    RESULT=np.zeros(X.shape)
#    print(X)
    RESULT[:-1,:-1]=((s[:-1,1:]-s[:-1,:-1])/(X[:-1,1:]-X[:-1,:-1]))*((m[1:,:-1]-m[:-1,:-1])/(Z[1:,:-1]-Z[:-1,:-1]))-((m[:-1,1:]-m[:-1,:-1])/(X[:-1,1:]-X[:-1,:-1]))*((s[1:,:-1]-s[:-1,:-1])/(Z[1:,:-1]-Z[:-1,:-1]))

    RESULT[:-1,:-1]=RESULT[:-1,:-1]/(np.sqrt(((m[1:,:-1]-m[:-1,:-1])/(Z[1:,:-1]-Z[:-1,:-1]))**2+((m[:-1,1:]-m[:-1,:-1])/(X[:-1,1:]-X[:-1,:-1]))**2)*np.sqrt(((s[1:,:-1]-s[:-1,:-1])/(Z[1:,:-1]-Z[:-1,:-1]))**2+((s[:-1,1:]-s[:-1,:-1])/(X[:-1,1:]-X[:-1,:-1]))**2))

    RESULT[:,-1]=RESULT[:,-2]
    RESULT[-1,:]=RESULT[-2,:]
    RESULT=RESULT**2
    return RESULT

#dir_comparison=sys.argv[1]
#for i in range(4):
#cross_section="track"+str(i)+".profile"
cross_section=sys.argv[1]
x1,z1,v1=np.loadtxt(cross_section,unpack=True)

#cross_section_with_constraint=dir_comparison+"/track"+str(i)+".profile"
cross_section_with_constraint=sys.argv[2]
x2,z2,v2=np.loadtxt(cross_section_with_constraint,unpack=True)
z=np.loadtxt("zpoints.txt",unpack=True)
nz=z.size
n_points=x1.size
nx=n_points//nz
print(f"nx={nx}, nz={nz}")
X=np.reshape(x1,(nz,nx),order='C')
Z=np.reshape(z1,(nz,nx),order='C')
V1=np.reshape(v1,(nz,nx),order='C')
V2=np.reshape(v2,(nz,nx),order='C')
#sim=similarity.compute_sim(V1,V2,win_size,method=comparison_method)
sim=cross_gradient(X,Z,V1,V2)
sim_flat=sim.flatten()
print(f"Max similarity={np.max(sim_flat)}, min similarity={np.min(sim_flat)}")
print("")
data=np.hstack((x1[:,np.newaxis],z1[:,np.newaxis],sim_flat[:,np.newaxis]))
#out_file="sim"+str(i)+".profile"
out_file=sys.argv[3]
np.savetxt(out_file,data,fmt="%25.6f")



