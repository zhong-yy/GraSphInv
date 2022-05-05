import sys
import numpy as np

lats, lons, residuals = np.loadtxt(sys.argv[1], unpack=True)
mean=np.mean(residuals)
std=np.std(residuals)
print("Mean=%.2g, STD=%.2g"%(mean,std))
maxv=np.max(residuals)
minv=np.min(residuals)
#print(maxv)
#print(minv)

