import sys
import numpy as np

x=np.zeros([2])
y=np.zeros([2])
x[0]=sys.argv[1]
x[1]=sys.argv[2]
y[0]=sys.argv[3]
y[1]=sys.argv[4]
f1=open(sys.argv[5],'a')#append

print(y)
print(x)

npoints=100
xs=np.concatenate( ( np.linspace(x[0],x[1],npoints),np.ones([npoints])*x[1],np.linspace(x[1],x[0],npoints),np.ones([npoints])*x[0]) )
ys=np.concatenate( ( np.ones([npoints])*y[0],np.linspace(y[0],y[1],npoints),np.ones([npoints])*y[1],np.linspace(y[1],y[0],npoints)) )
f1.write(">")
for i in np.arange(np.size(ys)):
	f1.write(" %.12e  %.12e\r\n" % (xs[i], ys[i]))
