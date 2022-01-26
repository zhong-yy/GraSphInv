import numpy as np
import scipy as sp
from numba import jit

def normalized_cross_correlation(a,b):
    ncc=(np.sum(a*b)+1e-5)/(np.sqrt(np.sum(a**2))*np.sqrt(np.sum(b**2))+1e-5)
    return ncc

def pearson_corr_coef(a,b):
    a2=a-np.mean(a)
    b2=b-np.mean(b)
    cc=np.sum(a2*b2)/(np.sqrt(np.sum(a2**2))*np.sqrt(np.sum(b2**2))+1e-9)
    return cc
def square_diff(a,b):
    result=np.sum((a-b)**2)/(np.sqrt(np.sum(a**2))*np.sqrt(np.sum(b**2))+1e-9)
    return result
def norm_diff(a,b):
    a2=(a-np.min(a))/(np.max(a)-np.min(a))
    b2=(b-np.min(b))/(np.max(b)-np.min(b))
    return (a2-b2)**2
    
#https://en.wikipedia.org/wiki/Structural_similarity
def ssim(a,b): #https://blog.csdn.net/hedgehog__/article/details/107257755
    mu_a=np.mean(a)
    mu_b=np.mean(b)
    std_a=np.std(a,ddof=1)
    std_b=np.std(b,ddof=1)
    std_ab=1/(np.size(a)-1)*np.sum((a-mu_a)*(b-mu_b))
    C1=1e-4
    C2=1e-4
    C3=C2/2
    I=(2*mu_a*mu_b+C1)/(mu_a*mu_a+mu_b*mu_b+C1)
    c=(2*std_a*std_b+C2)/(std_a*std_a+std_b*std_b+C2)
    s=(std_ab+C3)/(std_a*std_b+C3)
    #result=I*c*s
    result=(2*mu_a*mu_b+C1)*(2*std_ab+C2)/((mu_a*mu_a+mu_b*mu_b+C1)*(std_a*std_a+std_b*std_b+C2))
    return result


#sliding templ through image
def matchTemp(image,templ,method=normalized_cross_correlation):
    ny_im,nx_im=image.shape
    ny_te,nx_te=templ.shape
    result=np.zeros(image.shape)
    dxL=nx_te//2
    dxR=nx_te-dxL-1
    dyB=ny_te//2
    dyT=ny_te-dyB-1 
    image_pad=np.zeros((ny_im+dyB+dyT,nx_im+dxL+dxR))
    image_pad[dyB:ny_im+dyB,dxL:nx_im+dxL]=image
    for i in range(ny_im):
        for j in range(nx_im):
            result[i,j]=0
            result[i,j]=method(image_pad[i:i+ny_te,j:j+nx_te],templ)
    return result
def match_and_compute_sim(a,b,window_size,method=normalized_cross_correlation):
    #assert a.shape==b.shape, "The two images cannot be compared because they have different sizes."
    ny,nx=a.shape
    height=window_size[0]
    width=window_size[1]    
    result=np.zeros(a.shape)
    for i in range(0,ny-height,60):
        for j in range(0,nx-width,60):
            match_result=matchTemp(b,a[i:i+height,j:j+width],method=method)
            result[i:i+height,j:j+width]=np.max(match_result)
    return result
def compute_sim(a,b,window_size,method=normalized_cross_correlation):
    assert a.shape==b.shape, "The two images cannot be compared because they have different sizes."
    ny,nx=a.shape
    height=window_size[0]
    width=window_size[1]
    dxL=width//2
    dxR=width-dxL-1
    dyB=height//2
    dyT=height-dyB-1
    a_pad=np.zeros((ny+dyB+dyT,nx+dxL+dxR))
    b_pad=np.zeros((ny+dyB+dyT,nx+dxL+dxR))
    a_pad[dyB:ny+dyB,dxL:nx+dxL]=a
    b_pad[dyB:ny+dyB,dxL:nx+dxL]=b  
    result=np.zeros(a.shape)
    for i in range(0,ny):
        for j in range(0,nx):
            result[i,j]=method(a_pad[i:i+height,j:j+width],b_pad[i:i+height,j:j+width])
    print("Average similarity ",np.mean(result))
    return result

