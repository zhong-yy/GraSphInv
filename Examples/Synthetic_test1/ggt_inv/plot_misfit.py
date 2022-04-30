import numpy as np
import matplotlib.pyplot as plt
import string 
def plot_lambda_misfit(filename_iteration,filename_lambda,outname):
    data=np.loadtxt(filename_iteration,skiprows=1)
    data2=np.loadtxt(filename_lambda,skiprows=1)
    iteration=data[:,0]
    num_lam=data[:,1]
    lam=data2[:,0]
    misfit=data2[:,1]

    iteration_number=iteration[-1]
    print(iteration_number)
    plt.figure(figsize=[4,3]) #单位inch
    start=0
    end=int(num_lam[0])
    mk=('v','P','o')
    for i in range(0,int(iteration_number)):
        #print(ls[i])
        end=int(start+num_lam[i])    
        plt.loglog(lam[start:end],misfit[start:end],label='Iteration'+str(i+1),color='black',linestyle='-',marker=mk[i],markersize=5)
        plt.xlabel('Regularization parameter '+r'$\lambda$',fontsize=10)
        plt.ylabel('Misfit',fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.grid(b=True,axis='y',which='both')
        
        plt.legend(frameon=True,fancybox=True)
        print([start,end])        
        start=int(start+num_lam[i])
    plt.gca().invert_xaxis()

        
    plt.savefig(outname,dpi=300,bbox_inches='tight')
    
def plot_misfit(filename,outname):
    data=np.loadtxt(filename,skiprows=1)
    iteration=data[:,0]
    misfit=data[:,2]
    
    plt.figure(figsize=[6,5])
    plt.semilogy(iteration,misfit,'.-',label='L-curve')
    plt.xlabel('Iteration')
    plt.ylabel('misfit')

    plt.savefig(outname,dpi=300,bbox_inches='tight')
    
plot_lambda_misfit('Iteration_misfit_GN','lambda_misfit_GN','lambda_misfit.jpg')
plot_misfit('Iteration_misfit_GN','misfit.jpg')


