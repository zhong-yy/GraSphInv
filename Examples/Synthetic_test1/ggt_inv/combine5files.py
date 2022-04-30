import sys
import numpy as np

#G=np.loadtxt('dobs_g_r')
TRR=np.loadtxt('dobs_T_rr')
TRTHETA=np.loadtxt('dobs_T_rtheta')
TRPHI=np.loadtxt('dobs_T_rphi')
TTHETATHETA=np.loadtxt('dobs_T_thetatheta')
TTHETAPHI=np.loadtxt('dobs_T_thetaphi')
#TPHIPHI=np.loadtxt('dobs_T_phiphi')

rows,cols=TRR.shape

NEW=np.zeros((rows,7))

#latitude and longitude
NEW[:,[0,1]]=TRR[:,[0,1]]

#Trr
NEW[:,2]=TRR[:,2]

#Trtheta
NEW[:,3]=TRTHETA[:,2]
NEW[:,4]=TRPHI[:,2]
NEW[:,5]=TTHETATHETA[:,2]
NEW[:,6]=TTHETAPHI[:,2]

np.savetxt('dobs_FTG',NEW,'%12.7f %12.7f %23.15e %23.15e %23.15e %23.15e %23.15e')



