import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

fts=9
cm = 1/2.54
fig,ax=plt.subplots(figsize=(7.5*cm,7.5*cm))
lon,lat,r,v=np.loadtxt("./intp_vol_data_ada.txt",usecols=[6,7,8,10],unpack=True)
lon,lat,r,v_non_ada=np.loadtxt("./intp_vol_data_non_ada.txt",usecols=[6,7,8,10],unpack=True)
x=np.linspace(-155,155,2)
#plt.scatter(v,v_non_ada,c=6379.137-r)
ax.plot(v,v_non_ada,'o',markeredgecolor='black',markerfacecolor='white',markeredgewidth=0.3,markersize=2.5)
#ax.scatter(v,v_non_ada)
#ax.scatter(v,v_non_ada,marker='o',markeredgecolor='black',markerfacecolor='white',fillstyle='none',markersize=2)
ax.plot(x,x,color='blue')
ax.set_xlabel('Inverted density values from\n the adaptively refined mesh (kg/m$^3$)',fontsize=fts)
ax.set_ylabel('Inverted density values from\n the uniform mesh (kg/m$^3$)',fontsize=fts)
ax.set_aspect('equal')
ax.tick_params(which='both',labelsize=fts)
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(10))
#ax.set_xticks(np.arange(-100,150,50))
#ax.set_yticks(np.arange(-100,150,50))
ax.set_xlim([-155,155])
ax.set_ylim([-155,155])
ax.grid()
plt.savefig('cross_plot_2.jpg',dpi=300,bbox_inches='tight')
