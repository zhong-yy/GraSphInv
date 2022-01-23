import sys
import numpy as np

#load data
#https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html

#col1(latitude) col2(longitude) col3(vertical gravity with positive axis upward)
G_bouguer= np.loadtxt("Bouguer_g_r")

G_crust=np.loadtxt("../Crust_Correction/crystalline_crust_g_r")

G_sediments=np.loadtxt("../Crust_Correction/sediments_g_r")

G_moho=np.loadtxt("../Crust_Correction/moho_g_r")


G_new=G_bouguer.copy()
G_new[:,2]=G_new[:,2]-G_crust[:,2]-G_sediments[:,2]-G_moho[:,2]

#Usage of numpy.savetxt, please see:#https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
np.savetxt("mantle_g_r_EGM2008",G_new,'%20.7f %20.7f %30.15e')
