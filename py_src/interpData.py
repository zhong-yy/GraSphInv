import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import RBFInterpolator
from scipy.interpolate import LinearNDInterpolator
def geo2cart(lon,lat,r):
    lon=lon*np.pi/180.0
    lat=lat*np.pi/180.0
    x=r*np.cos(lat)*np.cos(lon)
    y=r*np.cos(lat)*np.sin(lon)
    z=r*np.sin(lat)
    return x,y,z
def interp_lon_slice2(filename,lon0,start_lat,stop_lat,num_lat,start_r,stop_r,num_r,cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
    # kernel: 'thin_plate_spline', 'linear', 'cubic'
    #Load data
    lon,lat,r,v=np.loadtxt(filename,skiprows=skiprows,usecols=cols,unpack=True)
    x_cart,y_cart,z_cart=geo2cart(lon,lat,r)
    min_x=np.min(x_cart)
    min_y=np.min(y_cart)
    min_z=np.min(z_cart)
    
    print(min_x,min_y,min_z)
    print(np.max(x_cart),np.max(y_cart),np.max(z_cart))
    
    x_cart=x_cart-min_x
    y_cart=y_cart-min_y
    z_cart=z_cart-min_z

    xyz=np.hstack((x_cart[:,np.newaxis],y_cart[:,np.newaxis],z_cart[:,np.newaxis]))

    #Construct interpolation grid points
    loni=np.array([lon0]) #lon
    lati=np.linspace(start_lat,stop_lat,num_lat) #lat
    ri=np.linspace(start_r,stop_r,num_r) #r
    LONI,LATI,RI=np.meshgrid(loni,lati,ri) #geographic coordinates
   
    # cartesian coordinates
    XI,YI,ZI=geo2cart(LONI,LATI,RI)
    XI=XI-min_x
    YI=YI-min_y
    ZI=ZI-min_z #coordinate transformation is performed to stabilized the interpolation

    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()

    #Perform radial basis function interpolation
    #rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    #XYZI_flat=XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    #VI_flat=rbf(XYZI_flat)
    #VI=np.reshape(VI_flat,XI.shape)
    VI=griddata(xyz,v,(XI,YI,ZI),method='linear')

    LATI_2d=LATI[:,0,:] # shape of YI is (ny,nx,nz)
    RI_2d=RI[:,0,:] # shape of ZI is (ny,nx,nz)
    VI_2d=VI[:,0,:]
    return LATI_2d,RI_2d,VI_2d

def interp_lon_slice(filename,lon0,start_lat,stop_lat,num_lat,start_r,stop_r,num_r,cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
    #Load data
    lon,lat,r,v=np.loadtxt(filename,skiprows=skiprows,usecols=cols,unpack=True)
    x_cart,y_cart,z_cart=geo2cart(lon,lat,r)
    min_x=np.min(x_cart)
    min_y=np.min(y_cart)
    min_z=np.min(z_cart)
    
    print(min_x,min_y,min_z)
    print(np.max(x_cart),np.max(y_cart),np.max(z_cart))
    
    x_cart=x_cart-min_x
    y_cart=y_cart-min_y
    z_cart=z_cart-min_z

    xyz=np.hstack((x_cart[:,np.newaxis],y_cart[:,np.newaxis],z_cart[:,np.newaxis]))

    #Construct interpolation grid points
    loni=np.array([lon0]) #lon
    lati=np.linspace(start_lat,stop_lat,num_lat) #lat
    ri=np.linspace(start_r,stop_r,num_r) #r
    LONI,LATI,RI=np.meshgrid(loni,lati,ri) #geographic coordinates
   
    # cartesian coordinates
    XI,YI,ZI=geo2cart(LONI,LATI,RI)
    XI=XI-min_x
    YI=YI-min_y
    ZI=ZI-min_z #coordinate transformation is performed to stabilized the interpolation

    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()

    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)

    LATI_2d=LATI[:,0,:] # shape of YI is (ny,nx,nz)
    RI_2d=RI[:,0,:] # shape of ZI is (ny,nx,nz)
    VI_2d=VI[:,0,:]
    return LATI_2d,RI_2d,VI_2d

def interp_lat_slice(filename,lat0,start_lon,stop_lon,num_lon,start_r,stop_r,num_r,cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
    #Load data
    lon,lat,r,v=np.loadtxt(filename,skiprows=skiprows,usecols=cols,unpack=True)
    x_cart,y_cart,z_cart=geo2cart(lon,lat,r)
    min_x=np.min(x_cart)
    min_y=np.min(y_cart)
    min_z=np.min(z_cart)
    
    print(min_x,min_y,min_z)
    print(np.max(x_cart),np.max(y_cart),np.max(z_cart))
    
    x_cart=x_cart-min_x
    y_cart=y_cart-min_y
    z_cart=z_cart-min_z

    xyz=np.hstack((x_cart[:,np.newaxis],y_cart[:,np.newaxis],z_cart[:,np.newaxis]))

    #Construct interpolation grid points
    lati=np.array([lat0])
    loni=np.linspace(start_lon,stop_lon,num_lon)
    ri=np.linspace(start_r,stop_r,num_r)
    LONI,LATI,RI=np.meshgrid(loni,lati,ri) #geographic coordinates
   
    # cartesian coordinates
    XI,YI,ZI=geo2cart(LONI,LATI,RI)
    XI=XI-min_x
    YI=YI-min_y
    ZI=ZI-min_z #coordinate transformation is performed to stabilized the interpolation

    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()

    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)

    LONI_2d=LONI[0,:,:] # shape of XI is (ny,nx,nz)
    RI_2d=RI[0,:,:] # shape of ZI is (ny,nx,nz)
    VI_2d=VI[0,:,:]
    return LONI_2d,RI_2d,VI_2d

def interp_r_slice(filename,r0,start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
    #Load data
    lon,lat,r,v=np.loadtxt(filename,skiprows=skiprows,usecols=cols,unpack=True)
    x_cart,y_cart,z_cart=geo2cart(lon,lat,r)
    min_x=np.min(x_cart)
    min_y=np.min(y_cart)
    min_z=np.min(z_cart)
    
    print(min_x,min_y,min_z)
    print(np.max(x_cart),np.max(y_cart),np.max(z_cart))
    
    x_cart=x_cart-min_x
    y_cart=y_cart-min_y
    z_cart=z_cart-min_z

    xyz=np.hstack((x_cart[:,np.newaxis],y_cart[:,np.newaxis],z_cart[:,np.newaxis]))

    #Construct interpolation grid points
    ri=np.array([r0])
    loni=np.linspace(start_lon,stop_lon,num_lon)
    lati=np.linspace(start_lat,stop_lat,num_lat)
    LONI,LATI,RI=np.meshgrid(loni,lati,ri) #geographic coordinates
   
    # cartesian coordinates
    XI,YI,ZI=geo2cart(LONI,LATI,RI)
    XI=XI-min_x
    YI=YI-min_y
    ZI=ZI-min_z #coordinate transformation is performed to stabilized the interpolation

    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()

    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)

    LONI_2d=LONI[:,:,0] # shape of XI is (ny,nx,nz)
    LATI_2d=LATI[:,:,0] # shape of ZI is (ny,nx,nz)
    VI_2d=VI[:,:,0]
    return LONI_2d,LATI_2d,VI_2d

def interp_r_slice2(filename,r0,start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
    #Load data
    lon,lat,r,v=np.loadtxt(filename,skiprows=skiprows,usecols=cols,unpack=True)
    x_cart,y_cart,z_cart=geo2cart(lon,lat,r)
    min_x=np.min(x_cart)
    min_y=np.min(y_cart)
    min_z=np.min(z_cart)
    
    print(min_x,min_y,min_z)
    print(np.max(x_cart),np.max(y_cart),np.max(z_cart))
    
    x_cart=x_cart-min_x
    y_cart=y_cart-min_y
    z_cart=z_cart-min_z

    xyz=np.hstack((x_cart[:,np.newaxis],y_cart[:,np.newaxis],z_cart[:,np.newaxis]))

    #Construct interpolation grid points
    ri=np.array([r0])
    loni=np.linspace(start_lon,stop_lon,num_lon)
    lati=np.linspace(start_lat,stop_lat,num_lat)
    LONI,LATI,RI=np.meshgrid(loni,lati,ri) #geographic coordinates
   
    # cartesian coordinates
    XI,YI,ZI=geo2cart(LONI,LATI,RI)
    XI=XI-min_x
    YI=YI-min_y
    ZI=ZI-min_z #coordinate transformation is performed to stabilized the interpolation

    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()

    #Perform radial basis function interpolation
    VI=griddata(xyz,v,(XI,YI,ZI))

    LONI_2d=LONI[:,:,0] # shape of XI is (ny,nx,nz)
    LATI_2d=LATI[:,:,0] # shape of ZI is (ny,nx,nz)
    VI_2d=VI[:,:,0]
    return LONI_2d,LATI_2d,VI_2d

def interp_volume(filename,start_lon,stop_lon,num_lon,start_lat,stop_lat,num_lat,start_r,stop_r,num_r,cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
    #Load data
    lon,lat,r,v=np.loadtxt(filename,skiprows=skiprows,usecols=cols,unpack=True)
    x_cart,y_cart,z_cart=geo2cart(lon,lat,r)
    min_x=np.min(x_cart)
    min_y=np.min(y_cart)
    min_z=np.min(z_cart)
    
    print(min_x,min_y,min_z)
    print(np.max(x_cart),np.max(y_cart),np.max(z_cart))
    print("max v=%f, min v=%f"%(np.max(v),np.min(v))) 
    
    x_cart=x_cart-min_x
    y_cart=y_cart-min_y
    z_cart=z_cart-min_z

    xyz=np.hstack((x_cart[:,np.newaxis],y_cart[:,np.newaxis],z_cart[:,np.newaxis]))

    #Construct interpolation grid points
    ri=np.linspace(start_r,stop_r,num_r) #r
    loni=np.linspace(start_lon,stop_lon,num_lon)
    lati=np.linspace(start_lat,stop_lat,num_lat)
    LONI,LATI,RI=np.meshgrid(loni,lati,ri) #geographic coordinates
   
    # cartesian coordinates
    XI,YI,ZI=geo2cart(LONI,LATI,RI)
    XI=XI-min_x
    YI=YI-min_y
    ZI=ZI-min_z #coordinate transformation is performed to stabilized the interpolation

    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()

    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)
    print("max v=%f, min v=%f"%(np.max(VI),np.min(VI))) 

    return LONI,LATI,RI,VI
    
def interp_volume_and_write_out(filename,
                                output_filename,
                                interp_dlon,interp_dlat,interp_dr,
                                start_lon,stop_lon,num_lon,
                                start_lat,stop_lat,num_lat,
                                start_r,stop_r,num_r,
                                cols=[6,7,8,10],skiprows=2,neighbors=None,kernel='linear'):
    #Load data
    lon,lat,r,v=np.loadtxt(filename,skiprows=skiprows,usecols=cols,unpack=True)
    x_cart,y_cart,z_cart=geo2cart(lon,lat,r)
    min_x=np.min(x_cart)
    min_y=np.min(y_cart)
    min_z=np.min(z_cart)
    
    print(min_x,min_y,min_z)
    print(np.max(x_cart),np.max(y_cart),np.max(z_cart))
    print("max v=%f, min v=%f"%(np.max(v),np.min(v))) 
    
    x_cart=x_cart-min_x
    y_cart=y_cart-min_y
    z_cart=z_cart-min_z

    xyz=np.hstack((x_cart[:,np.newaxis],y_cart[:,np.newaxis],z_cart[:,np.newaxis]))

    #Construct interpolation grid points
    ri=np.linspace(start_r,stop_r,num_r) #r
    loni=np.linspace(start_lon,stop_lon,num_lon)
    lati=np.linspace(start_lat,stop_lat,num_lat)
    LONI,LATI,RI=np.meshgrid(loni,lati,ri) #geographic coordinates
   
    # cartesian coordinates
    XI,YI,ZI=geo2cart(LONI,LATI,RI)
    XI=XI-min_x
    YI=YI-min_y
    ZI=ZI-min_z #coordinate transformation is performed to stabilized the interpolation

    XI_flat=XI.flatten()
    YI_flat=YI.flatten()
    ZI_flat=ZI.flatten()

    #Perform radial basis function interpolation
    rbf=RBFInterpolator(xyz,v,neighbors=neighbors,kernel=kernel)
    XYZI_flat=np.hstack((XI_flat[:,np.newaxis],YI_flat[:,np.newaxis],ZI_flat[:,np.newaxis]))
    VI_flat=rbf(XYZI_flat)
    VI=np.reshape(VI_flat,XI.shape)
    print("max v=%f, min v=%f"%(np.max(VI),np.min(VI)))
    
    lon_flat=LONI.flatten()
    lon_flat0=lon_flat-0.5*interp_dlon
    lon_flat1=lon_flat+0.5*interp_dlon    

    lat_flat=LATI.flatten()
    lat_flat0=lat_flat-0.5*interp_dlat            
    lat_flat1=lat_flat+0.5*interp_dlat            
    
    r_flat=RI.flatten()
    r_flat0=r_flat-0.5*interp_dr
    r_flat1=r_flat+0.5*interp_dr
    v_flat=VI.flatten()

    dep=6378.137-r_flat

    header="The first 6 parameters give the dimensions of a cell. The 7-9 parameters are the cell center. The 10th parameter is the value within the cell.\n%-22s%-24s%-24s%-24s%-24s%-24s%-24s%-24s%-24s%-24s%-24s"%("Longitude 0(degree)","Longitude 1(degree)","Latitude 0(degree)","Latitude 1(degree)","R 0(km)","R 1(km)","Longitude c(degree)","Latitude c(degree)","R c(km)","Depth (km)","value")
    #+"#Longitude 0(degree)   Longitude 1(degree)    #Latitude 0(degree)    Latitude 1(degree)     R 0(km)               R 1(km)                Longitude c(degree)    Latitude c(degree)     R c(km)                Depth (km)             value                  "
    data=np.hstack((lon_flat0[:,np.newaxis],lon_flat1[:,np.newaxis],lat_flat0[:,np.newaxis],lat_flat1[:,np.newaxis],r_flat0[:,np.newaxis],r_flat1[:,np.newaxis],lon_flat[:,np.newaxis],lat_flat[:,np.newaxis],r_flat[:,np.newaxis],dep[:,np.newaxis],v_flat[:,np.newaxis]))
    np.savetxt(output_filename,data,header=header,fmt="%-23.15e")

    return LONI,LATI,RI,VI

