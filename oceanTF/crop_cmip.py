import xarray as xr
from pyproj import Transformer
import os
import netCDF4 as nc
import numpy as np

def crop_cmip(cmipTfile,cmipSfile,roi,n_chunk=120):

    # output filenames
    outfileT = cmipTfile.replace('.nc','_cropped.nc')
    if os.path.exists(outfileT):
        os.remove(outfileT)
    outfileS = cmipSfile.replace('.nc','_cropped.nc')
    if os.path.exists(outfileS):
        os.remove(outfileS)
    
    # open datasets
    dsT = xr.open_dataset(cmipTfile)
    dsS = xr.open_dataset(cmipSfile)

    # get spatial coordinates
    lat = dsT['lat'].values
    lon = dsT['lon'].values
    latlon2polarstereo = Transformer.from_crs("EPSG:4326", "EPSG:3413", always_xy=True)
    x, y = latlon2polarstereo.transform(lon, lat)

    # get the indices for the region of interest
    x = x.reshape(x.shape[0]*x.shape[1])
    y = y.reshape(y.shape[0]*y.shape[1])
    roi_inds = (x > roi[0]) & (x < roi[1]) & (y > roi[2]) & (y < roi[3])

    # crop to roi
    x_roi = x[roi_inds]
    y_roi = y[roi_inds]
    lat_roi = lat.reshape(lat.shape[0]*lat.shape[1])[roi_inds]
    lon_roi = lon.reshape(lon.shape[0]*lon.shape[1])[roi_inds]

    # load depth and time coordinate
    lev = dsT['lev'].values
    t = dsT['time'].values

    # create the output datasets
    # temperature
    data_placeholder = np.full((len(t),len(lev),len(x_roi)), np.nan, dtype=np.float32)
    dsT_out = xr.Dataset(data_vars=dict(thetao=(["time", "lev", "points"], data_placeholder)),
                        coords=dict(time=t,lev=lev,x=("points", x_roi),y=("points", y_roi),lat=("points", lat_roi),lon=("points", lon_roi)),
                        attrs={"description": "Data cropped to roi, x-y coordinates in EPSG:3413"})
    dsT_out.to_netcdf(outfileT, encoding={"thetao": {"zlib": True, "complevel": 9}})
    # salinity
    dsS_out = xr.Dataset(data_vars=dict(so=(["time", "lev", "points"], data_placeholder)),
                        coords=dict(time=t,lev=lev,x=("points", x_roi),y=("points", y_roi),lat=("points", lat_roi),lon=("points", lon_roi)),
                        attrs={"description": "Data cropped to roi, x-y coordinates in EPSG:3413"})
    dsS_out.to_netcdf(outfileS, encoding={"so": {"zlib": True, "complevel": 9}})    

    # write the actual temperature data in chunks
    with nc.Dataset(outfileT, 'a') as nc_out:

        for start in range(0, len(t), n_chunk):
            end = min(start + n_chunk, len(t))
            print(f"Writing temperature time slice from {t[start]} to {t[end-1]}")

            # load data chunk and crop to roi
            thetao = dsT['thetao'][start:end,:,:,:].values
            thetao_roi = thetao.reshape(thetao.shape[0], thetao.shape[1], thetao.shape[2]*thetao.shape[3])[:,:,roi_inds]

            # write to netcdf
            nc_out.variables['thetao'][start:end, :, :] = thetao_roi

    # write the actual salinity data in chunks
    with nc.Dataset(outfileS, 'a') as nc_out:

        for start in range(0, len(t), n_chunk):
            end = min(start + n_chunk, len(t))
            print(f"Writing salinity time slice from {t[start]} to {t[end-1]}")

            # load data chunk and crop to roi
            so = dsS['so'][start:end,:,:,:].values
            so_roi = so.reshape(so.shape[0], so.shape[1], so.shape[2]*so.shape[3])[:,:,roi_inds]

            # write to netcdf
            nc_out.variables['so'][start:end, :, :] = so_roi

    # confirm finished
    print("All time slices written")