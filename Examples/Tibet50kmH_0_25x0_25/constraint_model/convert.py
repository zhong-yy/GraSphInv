import numpy as np
DATA=np.loadtxt('Tibetan_VsAbs')
DATA[:,3]=(DATA[:,3]-4500)/4500*100
np.savetxt('Tibetan_dVs_recomputed',DATA,'%20.7f %20.7f %20.7f %20.7f')
def relation(dvs,z):
    if dvs<=6 and dvs>0:
        drho=dvs*(7.3-z/100.0+dvs/4.0)
        if drho<0:
            print("drho=%f,dvs=%f,z=%f"%(drho,dvs,z))        
    elif dvs>6:
        drho=dvs*(8.8-z/100.0-7.0*(dvs-6.0)/40.0)
        if drho<0:
            print("drho=%f,dvs=%f,z=%f"%(drho,dvs,z))            
    else:
        drho=0
    return drho
DATA2=DATA.copy()
rows,cols=DATA.shape
for i in range(0,rows):
    DATA2[i,3]=relation(DATA[i,3],DATA[i,2])
np.savetxt('Tibetan_converted_rho',DATA2,'%20.7f %20.7f %20.7f %20.7f')


X=np.loadtxt('Tibetan_converted_rho')
del_rho=X[:,3]
print(X.shape)
print(np.max(del_rho))
print(np.min(del_rho))
print(np.max(DATA[:,3]))
print(np.min(DATA[:,3]))

