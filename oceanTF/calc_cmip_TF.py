import xarray as xr
import os
import numpy as np

def calc_cmip_TF(cmipTfile,cmipSfile,l1,l2,l3):

    # output file
    outfile = cmipTfile.replace('thetao','TF')
    if os.path.exists(outfile):
        os.remove(outfile)

    # open datasets
    dsT = xr.open_dataset(cmipTfile)
    dsS = xr.open_dataset(cmipSfile)

    # calculate TF
    z = dsT['lev'].broadcast_like(dsS['so'])/1e2 # nb convert from cm to m
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