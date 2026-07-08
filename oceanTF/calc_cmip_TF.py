import xarray as xr
import os
import numpy as np

def calc_cmip_TF(cmipTfile,cmipSfile,l1,l2,l3, depth_units='cm'):

    # output file
    outfile = cmipTfile.replace('thetao','TF')
    if os.path.exists(outfile):
        os.remove(outfile)

    # open datasets
    dsT = xr.open_dataset(cmipTfile)
    dsS = xr.open_dataset(cmipSfile)

    # calculate TF
    if depth_units=='cm':
        print('Depth units provided are interpreted as cm. Converting to m for calculation.')
        z = dsT['lev'].broadcast_like(dsS['so'])/1e2 # nb convert from cm to m
    elif depth_units=='m':
        print('Depth units provided are interpreted as m. No conversion applied.')
        z = dsT['lev'].broadcast_like(dsS['so'])
    else:
        print('Unrecognized depth units provided.  Please specify cm or m, or update function.')
        pass ## this will break it at next step...could write to throw error instead if feeling fancy
    Tfreeze = l1*dsS['so'] + l2 + l3*z
    TF = dsT['thetao'] - Tfreeze

    # set any TF<0 to 0, but preserve NaNs
    TF = TF.where((TF >= 0) | np.isnan(TF), 0)

    # save to netcdf
    ds_out = dsT.copy()
    ds_out = ds_out.drop_vars('thetao')
    TF.name = 'TF'
    ds_out['TF'] = TF.astype('float32')
    ds_out.to_netcdf(outfile, encoding={"TF": {"zlib": True, "complevel": 9, "dtype": 'float32'}})