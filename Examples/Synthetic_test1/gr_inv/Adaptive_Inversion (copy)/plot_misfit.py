import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
import string
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
def plot_lambda_misfit(filename_iteration,filename_lambda,outname):
    fts=9
    data=np.loadtxt(filename_iteration,skiprows=1)
    data2=np.loadtxt(filename_lambda,skiprows=1)
    iteration=data[:,0]
    num_lam=data[:,1]
    lam=data2[:,0]
    misfit=data2[:,1]

    iteration_number=iteration[-1]
    print(iteration_number)
    cm = 1/2.54
    plt.figure(figsize=(7.5*cm,5.5*cm)) #单位inch
    start=0
    end=int(num_lam[0])
    mk=('v','P','o','s')
    for i in range(0,int(iteration_number)):
        #print(ls[i])
        end=int(start+num_lam[i])    
        plt.loglog(lam[start:end],misfit[start:end],label='Iteration'+str(i+1),color='black',linestyle='-',marker=mk[i],markersize=5)
        for axis in [plt.gca().xaxis, plt.gca().yaxis]:
            formatter = FuncFormatter(lambda y, _: '{:.1f}'.format(y))
            axis.set_major_formatter(formatter)
        plt.gca().yaxis.set_minor_formatter(formatter)
        #plt.gca().ticklabel_format(style='plain',useMathText=True)
        plt.xlabel('Regularization parameter '+r'$\lambda$',fontsize=fts)
        plt.ylabel(r'Misfit $\sqrt{\frac{\Phi_d}{N_d}}$',fontsize=fts)
        plt.xticks(fontsize=fts)
        plt.yticks(fontsize=fts)
        plt.grid(b=True,axis='y',which='both')
        
        plt.legend(frameon=False,fancybox=False,fontsize=fts)
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
    
plot_lambda_misfit('Iteration_misfit_GN','lambda_misfit_GN','lambda_misfit.eps')
plot_misfit('Iteration_misfit_GN','misfit.jpg')


