import sys
import numpy as np

#load data
#https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html

#col1(latitude) col2(longitude) col3(vertical gravity with positive axis upward)
G_obs= np.loadtxt(sys.argv[1])


#swap the 1st and the 2nd colomuns (latitude and longitude)
G_obs[:,[0,1]]=G_obs[:,[1,0]]

#now, the format is
#col1(longitude) col2(latitude) col3(vertical gravity with positive axis upward)

#change the positive axis downward
G_obs[:,2]=-G_obs[:,2]

#now, the format is
#col1(longitude) col2(latitude) col3(vertical gravity with positive axis downward)
suffix=sys.argv[2]

#Save processed data
#Usage of numpy.savetxt, please see:#https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
np.savetxt(sys.argv[1]+suffix,G_obs,'%20.7f %20.7f %30.15e')

