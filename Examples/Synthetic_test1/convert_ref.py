import numpy as np
DATA=np.loadtxt('crg_model')
#DATA[:,3]=(DATA[:,3]-4500)/4500*100
print(np.max(DATA[:,3]))
print(np.min(DATA[:,3]))
def relation(dvs,z):
    if dvs<=6 and dvs>0:
        drho=dvs*(7.3-z/100.0+dvs/4.0)
    elif dvs>6:
        drho=dvs*(8.8-z/100.0-7.0*(dvs-6.0)/40.0)
    else:
        drho=0
    return drho
DATA2=DATA.copy()
rows,cols=DATA.shape
for i in range(0,rows):
    DATA2[i,3]=relation(DATA[i,3],DATA[i,2])
np.savetxt('ref_model2',DATA2,'%20.7f %20.7f %20.7f %20.7f')



X=np.loadtxt('ref_model')
del_rho=X[:,3]
print(X.shape)
print(relation(15,190))
print(np.max(del_rho))
print(np.min(del_rho))

