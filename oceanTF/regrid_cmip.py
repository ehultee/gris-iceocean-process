import xarray as xr
import os
from scipy.interpolate import griddata
import numpy as np
from pyproj import Transformer
from tqdm import tqdm

def regrid_cmip(cmipTFfile_cropped,xreg,yreg,zreg, which_variable='TF'):
    '''
    which_variable: str, optional. Default 'TF'.
        Must match the name of a variable in the dataset you want to regrid.
    '''

    # output filename
    outfile = cmipTFfile_cropped.replace('.nc','_regrid.nc')
    if os.path.exists(outfile):
        os.remove(outfile)

    # read in dataset
    ds = xr.open_dataset(cmipTFfile_cropped)
    t = ds['time'].values
    x = ds['x'].values
    y = ds['y'].values
    z = ds['lev'].values/1e2 # convert from cm to m (specific to CESM2-WACCM?)
    TF = ds[which_variable].values ## based on a function for regridding TF, but variable defined when called
    
    # first interpolate cmip variable onto regular depth grid
    TF_zreg = np.full((TF.shape[0], len(zreg), TF.shape[2]), np.nan)
    for i in range(TF_zreg.shape[0]):
        for j in range(TF_zreg.shape[2]):
            TF_zreg[i,:,j] = np.interp(zreg, z, TF[i,:,j])

    # now interpolate cmip variable onto regular horizontal grid
    TF_zreg_xyreg = np.full((TF.shape[0], len(zreg), len(yreg), len(xreg)), np.nan)
    points = np.column_stack((x, y))
    Xreg, Yreg = np.meshgrid(xreg, yreg)
    points_reg = np.column_stack((Xreg.ravel(), Yreg.ravel()))
    for i in tqdm(range(TF_zreg.shape[0])):
        #print(i)
        for k in range(TF_zreg.shape[1]):
            TF_zreg_xyreg[i,k,:,:] = griddata(points,TF_zreg[i,k,:],points_reg,method='linear').reshape(len(yreg),len(xreg))

    # lat-lon of regular grid
    polarstereo2latlon = Transformer.from_crs("EPSG:3413","EPSG:4326",always_xy=True)
    lonreg,latreg = polarstereo2latlon.transform(Xreg,Yreg)

    # catch any TF<0
    if which_variable=='TF':
        TF_zreg_xyreg[TF_zreg_xyreg<0] = 0
    else:
        print('Regridding non-TF variable. Did not correct negative values.')

    # put into xarray dataset
    ds_out = xr.Dataset({which_variable: (['time', 'depth', 'y', 'x'], TF_zreg_xyreg.astype(np.float32))},
                        coords={'time': t,'depth': zreg.astype(np.float32),
                                'y': yreg.astype(np.float32),'x': xreg.astype(np.float32),
                                'lat': (('y', 'x'), latreg.astype(np.float32)),'lon': (('y', 'x'), lonreg.astype(np.float32))})

    # save out to netcdf
    ds_out.to_netcdf(outfile, encoding={which_variable: {"zlib": True, "complevel": 9}})