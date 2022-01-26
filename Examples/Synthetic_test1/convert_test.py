import numpy as np
import matplotlib.pyplot as plt
#DATA=np.loadtxt('Tibetan_VsAbs')
#DATA[:,3]=(DATA[:,3]-4500)/4500*100
#print(np.max(DATA[:,3]))
#print(np.min(DATA[:,3]))
#np.savetxt('Tibetan_dVs_recomputed',DATA,'%20.7f %20.7f %20.7f %20.7f')
def relation(dvs,z):
    if dvs<=6 and dvs>0:
        drho=dvs*(7.3-z/100.0+dvs/4.0)
    elif dvs>6:
        drho=dvs*(8.8-z/100.0-7.0*(dvs-6.0)/40.0)
    else:
        drho=0
    return drho
nx=100
nz=9
x=np.linspace(1,25,nx)
z=np.linspace(80,200,nz)
#print(x.size)
y=np.zeros((nz,nx))
for j in range(0,z.size):
    for i in range(0,x.size):
        y[j,i]=relation(x[i],z[j])
    plt.plot(x,y[j,:],'-',label="z="+str(z[j])+"km")
print(relation(15,120))
plt.grid()
plt.ylabel(r"$\Delta\rho$ (kg/m$^3$)")
plt.xlabel(r"$v_s$ perturbation (%)")
plt.legend()
plt.show()
#print(y)
#DATA2=DATA.copy()
#rows,cols=DATA.shape
#for i in range(0,rows):
#    DATA2[i,3]=relation(DATA[i,3],DATA[i,2])
#np.savetxt('Tibetan_converted_rho',DATA2,'%20.7f %20.7f %20.7f %20.7f')


#X=np.loadtxt('Tibetan_converted_rho')
#del_rho=X[:,3]
#print(X.shape)
#print(np.max(del_rho))
#print(np.min(del_rho))

