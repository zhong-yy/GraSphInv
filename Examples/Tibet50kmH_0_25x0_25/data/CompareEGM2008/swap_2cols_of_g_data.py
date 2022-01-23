import sys
import numpy as np

#load data
#https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html

#col1(latitude) col2(longitude) col3(vertical gravity with positive axis upward)
G_obs= np.loadtxt(sys.argv[1])
G_predicted = np.loadtxt(sys.argv[2])

#swap the 1st and the 2nd colomuns (latitude and longitude)
G_obs[:,[0,1]]=G_obs[:,[1,0]]
G_predicted[:,[0,1]]=G_predicted[:,[1,0]]
#now, the format is
#col1(longitude) col2(latitude) col3(vertical gravity with positive axis upward)

#change the positive axis downward
G_obs[:,2]=-G_obs[:,2]
G_predicted[:,2]=-G_predicted[:,2]
#now, the format is
#col1(longitude) col2(latitude) col3(vertical gravity with positive axis downward)

#calculate residual
G_residual=G_obs.copy()
G_residual[:,2]-=G_predicted[:,2]

suffix=sys.argv[3]
#Save processed data
#Usage of numpy.savetxt, please see:#https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
np.savetxt(sys.argv[1]+suffix,G_obs,'%20.7f %20.7f %30.15e')
np.savetxt(sys.argv[2]+suffix,G_predicted,'%20.7f %20.7f %30.15e')
np.savetxt('residual'+suffix,G_residual,'%20.7f %20.7f %30.15e')

