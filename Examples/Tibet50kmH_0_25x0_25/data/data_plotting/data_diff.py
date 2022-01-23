import sys
import numpy as np

#load data
#https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html

#col1(latitude) col2(longitude) col3(vertical gravity with positive axis upward)
data1= np.loadtxt(sys.argv[1])
data2= np.loadtxt(sys.argv[2])

differences=data1.copy()
differences[:,2]=data1[:,2]-data2[:,2]
print(data1.shape)
print(differences.shape)
#Save processed data
#Usage of numpy.savetxt, please see:#https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
np.savetxt(sys.argv[3],differences,'%20.7f %20.7f %30.15e')

