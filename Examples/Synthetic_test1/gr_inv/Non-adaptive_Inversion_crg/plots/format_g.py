import sys
import numpy

lats, lons, g = numpy.loadtxt(sys.argv[1], unpack=True)
lats, lons, g_predicted = numpy.loadtxt(sys.argv[2], unpack=True)

suffix=sys.argv[3]

f1=open(sys.argv[1]+sys.argv[3],'w')
f2=open(sys.argv[2]+sys.argv[3],'w')
f3=open("residual"+sys.argv[3],'w')

for i in range(len(lats)):
	f1.write(" %.12e  %.12e %.12e \r\n" % (lons[i], lats[i],-g[i]))
	f2.write(" %.12e  %.12e %.12e \r\n" % (lons[i], lats[i],-g_predicted[i]))
	f3.write(" %.12e  %.12e %.12e \r\n" % (lons[i], lats[i],-(g[i]-g_predicted[i])))
