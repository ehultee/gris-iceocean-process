#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions 
@author: vincent
"""

import os
import sys
import copy
import numpy as np
import netCDF4 as nc
import csv
import pandas as pd
import pyproj
# import pwlf
from scipy.stats import norm
from scipy.interpolate import interp1d
# import statsmodels.tsa as tsa
# from statsmodels.tsa.arima.model import ARIMA

cwd = os.getcwd()+'/'

transformer1 = pyproj.Transformer.from_crs('epsg:4326','epsg:3413') #transformer of lat-lon to x-y(3413)
transformer2 = pyproj.Transformer.from_crs('epsg:3413','epsg:4326') #transformer of x-y(3413) to lat-lon

def freezingPoint(sal,depth):
    '''Freezing point temperature of seawater from Cowton et al. (2015)'''
    depth = abs(depth) #make sure to use depth values increasing downwards
    lb1 = -5.73e-2 #[°C psu-1] Cowton et al. (2015) Table 1
    lb2 = 8.32e-2 #[°C] Cowton et al. (2015) Table 1
    lb3 = -7.61e-4 #[°C m-1] Cowton et al. (2015) Table 1 (negative -> more negative as depth increases downwards)
    tempfr = lb1*sal+lb2+lb3*depth #Cowton et al. (2015) Eq.7
    return(tempfr)

def calcDpAveraged(values,depth,dmin,dmax):
    '''
    Calculates depth-averaged values
    '''
    if np.all(depth<=0): #function must use depth positive downwards
        depth = -1*depth
    if depth[0]>0:
        depth  = np.append(0,depth)
        values = np.append(values[0],values) #approxiumation for value at surface

    # Limit domain to [dmin:dmax] range # 
    dp0  = depth[depth<abs(dmax)]
    vl0  = values[depth<abs(dmax)]
    dp1  = dp0[dp0>abs(dmin)]
    vl1  = vl0[dp0>abs(dmin)]
    # Interpolate values at dmin and dmax #
    if np.any(depth<=abs(dmin)):
        dp1  = np.append(dmin,dp1)
        vl1  = np.append(linearInterpolation(depth,values,[dmin]),vl1)
    if np.any(depth>=abs(dmax)):
        dp1 = np.append(dp1,dmax)
        vl1 = np.append(vl1,linearInterpolation(depth,values,[dmax]))

    deltaz      = dp1[1:]-dp1[0:-1]
    valuesav    = (vl1[1:]+vl1[0:-1])/2 #mean value over interval
    averagedval = np.sum(deltaz*valuesav)/np.sum(deltaz) #values integrated over [dmin:dmax]
    return(averagedval)


def qmProjCannon2015(timeObs,srObs,timeMod,srMod,HistDates,windowProj=np.nan,ConstraintMin=np.nan,ConstraintMax=np.nan,exclSigma=np.nan):
    '''

    Working(12Oct2022):
        Note that Cannon et al. separate the RCP8.5 proj period in chunks of
        30 years to avoid large differences between min(deltaproj) and max(deltaproj)
        Here, I use a running CDF over windowProj years in Projection period

    Computes Quantile Mapping following method of Cannon et al. (2015)
    Observations over period HistDates[0] to HistDates[1] constrain the CDFs
    chunkProj allows to compute deltaproj between quantiles over separate chunks of the Projection Perio
        note that Cannon et al. (2015) separate the RCP8.5 in chunks of 30 years
    ConstraintMin can enforce the calibration inverse CDF to go through (0,ConstraintMin)
    ConstraintMax can enforce the calibration inverse CDF to go through (1,ConstraintMax)
    exclSigma sets values in final sr that are beyond +/-exclSigma std devs to +/-exclSigma std devs
    Returns:
        time of Historic period and associated corrected series
        time of Projection period and associated corrected series
    '''

    valsObsfull = np.copy(srObs)
    valsModfull = np.copy(srMod)

    ### Define historical period for Observations ###
    indsObs   = np.logical_and(timeObs>=HistDates[0],timeObs<=HistDates[1])
    tmobs     = timeObs[indsObs]
    valsObs   = valsObsfull[indsObs]

    ### Sort and estimate Observation CDF ###
    sortObs   = np.sort(valsObs)
    isortObs  = np.argsort(valsObs)
    cdfObs    = np.linspace(1,len(tmobs),len(tmobs))/len(tmobs)

    ### Define historical and projection periods for Model ###
    indshist = np.logical_and(timeMod>=HistDates[0],timeMod<=HistDates[1])
    if(np.any(timeMod>HistDates[1])):
        indsproj = (timeMod>HistDates[1])
        i0Proj   = np.where
    elif(np.any(timeMod<HistDates[0])):
        indsproj = (timeMod<HistDates[0])
        typeProj = 'preproj'
    i0Proj    = np.where(indsproj)[0][0]
    i1Proj    = np.where(indsproj)[0][-1]
    tmmodhist = timeMod[indshist]

    ### Sort and estimate Model hist CDF ###
    sortModhist   = np.sort(srMod[indshist])
    isortModhist  = np.argsort(srMod[indshist])
    cdfModhist    = np.linspace(1,len(tmmodhist),len(tmmodhist))/len(tmmodhist)

    ### Linear interpolation functions of Mod hist CDF values ###
    #fVtoP_Modhist  = interp1d(sortModhist,cdfModhist,kind='linear',fill_value='extrapolate')
    fPtoV_Modhist  = interp1d(cdfModhist,sortModhist,kind='linear',fill_value='extrapolate')

    ### Enforce constraints for the interpolation function of Obs CDF values ###
    cdfObsComp     = copy.deepcopy(cdfObs)
    sortObsComp    = copy.deepcopy(sortObs)
    # Function will assign ConstraintMin to all CDF values below cdfObs that corresponds to ConstraintMin #
    if(np.isnan(ConstraintMin)==False):
        # Remove Obs values below ConstraintMin #
        vmin   = copy.deepcopy(ConstraintMin)
        i0keep = np.where(sortObsComp>=ConstraintMin)[0][0]
        if(i0keep==0):
            cdfValForMin  = 0
        elif(i0keep>0): #interpolate value of the CDF at ConstraintMin
            cdfValForMin  = cdfObsComp[i0keep-1]+(ConstraintMin-sortObsComp[i0keep-1])*(cdfObsComp[i0keep]-cdfObsComp[i0keep-1])/(sortObsComp[i0keep]-sortObsComp[i0keep-1])
    else:
        # Linearly extrapolate values at CDF==0 #
        i0keep       = 0 #keep all indices
        cdfValForMin = 0
        vmin         = sortObsComp[0]-cdfObsComp[0]*(sortObsComp[1]-sortObsComp[0])/(cdfObsComp[1]-cdfObsComp[0])
    cdfObsComp  = np.append(cdfValForMin,cdfObsComp[i0keep:])
    sortObsComp = np.append(vmin,sortObsComp[i0keep:])
    # Function will assign ConstraintMax to all CDF values above cdfObs that corresponds to ConstraintMax #
    if(np.isnan(ConstraintMax)==False):
        # Remove Obs values above ConstraintMax #
        vmax   = copy.deepcopy(ConstraintMax)
        i1keep = np.where(sortObsComp<=ConstraintMax)[0][-1]
        if(i1keep==len(sortObsComp)-1):
            cdfValForMax = 1.0
        elif(i1keep<len(sortObsComp)-1): #interpolate value of the CDF at ConstraintMax
            cdfValForMax = cdfObsComp[i1keep]+(ConstraintMax-sortObsComp[i1keep])*(cdfObsComp[i1keep+1]-cdfObsComp[i1keep])/(sortObsComp[i1keep+1]-sortObsComp[i1keep])
    else:
        # Linearly extrapolate values at CDF==1 #
        i1keep       = len(sortObsComp)-1 #keep all indices
        cdfValForMax = 1
        vmax         = sortObsComp[-1]+(1-cdfObsComp[-1])*(sortObsComp[-1]-sortObsComp[-2])/(cdfObsComp[-1]-cdfObsComp[-2])
    cdfObsComp  = np.append(cdfObsComp[0:i1keep+1],cdfValForMax)
    sortObsComp = np.append(sortObsComp[0:i1keep+1],vmax)
    ### Fix computed values corresponding to CDF values extending beyond cdfValForMin and cdfValForMax ###
    OutOfBounds = (vmin,vmax)

    ### Linear interpolation function ###
    fPtoV_Obs = interp1d(cdfObsComp,sortObsComp,kind='linear',bounds_error=False,fill_value=OutOfBounds)

    ### Method of Cannon et al. (2015) ###
    deltahist    = 0
    xhathist     = fPtoV_Obs(cdfModhist)
    bchist0      = xhathist+deltahist

    lsisorttfModhist  = isortModhist.tolist()
    bchist1           = np.array([bchist0[lsisorttfModhist.index(ii)] for ii in range(sum(indshist))])

    ##### Now deal with projection period #####
    iinbproj       = np.where(indsproj)[0] #index number of the projection indices
    tmmodproj      = timeMod[indsproj] #projection period array only
    nwdw           = len(timeMod[timeMod<timeMod[0]+windowProj]) #nb of indices in windowProj
    bcproj0        = np.zeros(len(tmmodproj))
    valscdfModproj = np.zeros(len(tmmodproj))
    for projii,ii in enumerate(iinbproj):
        iivalorg = srMod[ii] #original value at ii
        iitm     = timeMod[ii] #time at ii
        indswdwproj = abs(iitm-tmmodproj)<=windowProj/2 #indices of proj within window
        if(np.sum(indswdwproj)<nwdw):
            if(indswdwproj[0]):
                #go from tmmodproj[0] to tmmodproj+windowProj
                indswdwproj[0:nwdw] = True
            elif(indswdwproj[-1]):
                #go from tmmodproj-windowProj to tmmodproj[-1]
                indswdwproj[-1*nwdw:] = True
        indswdwii   = np.zeros(len(timeMod)).astype(bool)
        indswdwii[i0Proj+np.where(indswdwproj)[0]] = True
        srmodproj   = srMod[indswdwii] #select values
        wdwtmproj   = timeMod[indswdwii] #select time
        posinsort   = np.where(np.sort(srmodproj)==iivalorg)[0][0] #position of ii in sorted values of window
        valscdfModproj[projii]  = (1+posinsort)/len(wdwtmproj) #value of ii in CDF of window
        #print(posinsort,valscdfModproj[projii])
    deltaproj = srMod[iinbproj]-fPtoV_Modhist(valscdfModproj) #change in quantile wrt hist
    xhatproj  = fPtoV_Obs(valscdfModproj) #values in Obs CDF
    bcproj0   = xhatproj+deltaproj
    # Ensure deltaproj did not move values beyond constraints #
    if(np.isnan(ConstraintMin)==False):
        bcproj0  = np.maximum(bcproj0,ConstraintMin)
    if(np.isnan(ConstraintMax)==False):
        bcproj0  = np.minimum(bcproj0,ConstraintMax)
    bcproj1 = np.array(bcproj0) #convert to numpy array

    ### Constrain values beyond exclSigma standard deviations if enforced ###
    if(np.isnan(exclSigma)==False):
        lowsigval  = np.mean(bchist1)-exclSigma*np.std(bchist1)
        upsigval   = np.mean(bchist1)+exclSigma*np.std(bchist1)
        bchist     = np.maximum(np.minimum(bchist1,upsigval),lowsigval)
        fullsr1   = np.append(bchist1,bcproj1) #order does not matter here
        lowsigval = np.mean(fullsr1)-exclSigma*np.std(fullsr1)
        upsigval  = np.mean(fullsr1)+exclSigma*np.std(fullsr1)
        bcproj1   = np.maximum(np.minimum(bcproj1,upsigval),lowsigval)

    ### Enforce continuity between annual means of hist and proj ###
    indsinproj    = np.arange(0,len(bcproj1),1)
    if(tmmodhist[-1]<tmmodproj[0]):
        # Proj period after Hist period #
        iyrhist     = np.floor(tmmodhist)==np.floor(tmmodhist[-1])
        iyrproj     = np.floor(tmmodproj)==np.floor(tmmodproj[0])
        #gap         = bchist[-1]-bcproj2[0]
        gap         = np.mean(bchist1[iyrhist])-np.mean(bcproj1[iyrproj])
        adjterm     = gap*((len(bcproj1)-indsinproj)/len(bcproj1))
    elif(tmmodhist[0]>tmmodproj[-1]):
        # Proj period before Hist period #
        iyrhist     = np.floor(tmmodhist)==np.floor(tmmodhist[0])
        iyrproj     = np.floor(tmmodproj)==np.floor(tmmodproj[-1])
        #gap         = bchist[0]-bcproj2[-1]
        gap         = np.mean(bchist1[iyrhist])-np.mean(bcproj1[iyrproj])
        adjterm     = gap*indsinproj/len(bcproj1)
    # Make sure constraints are still enforced #
    if(np.isnan(ConstraintMin)==False):
        adjterm = np.maximum(adjterm,ConstraintMin-bcproj1)
    if(np.isnan(ConstraintMax)==False):
        adjterm = np.minimum(adjterm,ConstraintMax-bcproj1)
    bcproj      = bcproj1+adjterm

    # Return #
    return(tmmodhist,bchist,tmmodproj,bcproj)
   
def getBadECCO2arcticPoints():
    '''
    Manual removal of ECCO2arctic points showing unrealistic TF series
    '''
    LocRemoval = f'{cwd}/InputOutput/ecco2arctic_removeinds.csv'
    data       = np.genfromtxt(LocRemoval,dtype=int,delimiter=',')
    indy,indx  = data[:,0],data[:,1]
    return(indy,indx)

def getxycoords(listoflatlonpairs):
    '''
    Returns x,y coordinates for a batch of lat-lon pairs
    Each lat array and lon array should be passed as 2D
    '''
    npairs  = int(len(listoflatlonpairs)/2)
    outlist = []
    for pp in range(npairs):
        latspp,lonspp = listoflatlonpairs[2*pp],listoflatlonpairs[2*pp+1]
        xpp,ypp       = np.zeros_like(latspp),np.zeros_like(latspp)
        for ii in range(np.shape(latspp)[0]):
            for jj in range(np.shape(latspp)[1]):
                (xx,yy) = transformer1.transform(latspp[ii,jj],lonspp[ii,jj])
                xpp[ii,jj] = xx
                ypp[ii,jj] = yy
        outlist.append(xpp)
        outlist.append(ypp)
    return(outlist)

def findNearestECCO2arcticPoint(latcoords,loncoords,bathyminECCO=0,NumNbrs=1):
    '''
    Returns the NumNbrs closest grid points of the ECCO2arctic grid from a given lat,lon pair
    latcoords and loncoords can be a single scalar or an array of values
    Output is in format of a 3D array:
        row for [0] closest neighbor [index(y),index(x),distance]
        row for [1] closest neighbor [index(y),index(x),distance]
        ...
        row for [NumNbrs-1] closest neighbor [index(y),index(x),distance]
        Each 3D entry corresponds to a pair from latcoords, loncoords
    '''
    
    DirECCO2arctic = f'{cwd}/InputOutput/'
    ECCOloc        = DirECCO2arctic+'tfECCO2arctic_only1992.nc'

    ### Manual removal ###
    removey,removex   = getBadECCO2arcticPoints()

    ### Load ECCO2arctic ###
    ds            = nc.Dataset(ECCOloc)
    depth         = np.array(ds.variables['depth'])
    lat           = np.array(ds.variables['lat'])
    lon           = np.array(ds.variables['lon'])
    tf            = np.array(ds.variables['thermalforcing'][0:2,:,:,:])
    ds.close()
    
    ### Remove bad ECCO2arctic points ###
    tf[:,:,removey,removex] = 1.1e20

    ### Coordinates in x,y,z ###
    cECCO = getxycoords([lat,lon])
    xecco,yecco = cECCO[0],cECCO[1]
    izz   = np.where(abs(depth)>=abs(bathyminECCO))[0][0]

    ### Convert scalar to single-element list if needed ###
    if(np.size(latcoords)==1):
        latcoords = [latcoords]
        loncoords = [loncoords]
      
    NumPoints = len(latcoords) #number of points for which we are looking for neighbors
    output    = np.zeros((NumPoints,NumNbrs,3))
    ### Find nearest neighbor for each input pair ###
    for kk in range(NumPoints):
        (xx,yy) = transformer1.transform(latcoords[kk],loncoords[kk])
        dists                       = np.sqrt((xx-xecco)**2+(yy-yecco)**2)
        dists[tf[0,izz,:,:]>1e10]   = 1.1e20
        ### Find indices of sorted distances ###
        isrt = np.unravel_index(np.argsort(dists,axis=None),dists.shape)
        for jj in range(NumNbrs):
            ### Find jjth smallest distance and corresponding indy,indx ###
            dmin   = dists[isrt[0][jj],isrt[1][jj]]
            indymin   = isrt[0][jj]
            indxmin   = isrt[1][jj]
            ### Save in output ###
            output[kk,jj,:] = np.array([indymin,indxmin,dmin])
    return(output)

def getMonthIndex(timearray):
    '''
    Gets the monthly index for each time step
    '''
    daysmonth        = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    daysmonthcum     = np.cumsum(daysmonth)
    daysmidmonthcum  = daysmonthcum-daysmonth/2
    midmonthfracyear = daysmidmonthcum/np.sum(daysmonth)
    fracyear         = timearray-np.floor(timearray)
    indsmonth        = np.array([np.argmin(abs(frac-midmonthfracyear)) for frac in fracyear])
    return(indsmonth)

def getMonthlyEffects(timearray,valsarray):
    '''
    Gets the mean monthly deviation from the annual mean of a time series
    '''
    ### Monthly indices ###
    indsmonth = getMonthIndex(timearray)
    ### Yearly averages ###
    yearlyMn  = np.zeros(len(timearray))
    yrs       = np.unique(np.floor(timearray))
    nyr       = len(yrs)
    for yy in range(nyr):
        yrinds           = np.where(np.floor(timearray)==yrs[yy])[0]
        yearlyMn[yrinds] = np.mean(valsarray[yrinds])

    ### Calculate monthly effects ###
    fullmontheffect  = valsarray-yearlyMn
    monthsignal      = np.zeros(12)
    for mm in range(12):
        inds            = np.where(indsmonth==mm)[0]
        monthsignal[mm] = np.mean(fullmontheffect[inds])
    ### Set full array of averaged monthly effects ###
    fullavgmonthsignal = np.array([monthsignal[indsmonth[jj]] for jj in range(len(timearray))])
    return(fullavgmonthsignal)

def getDevFromYearlyMean(timearray,valsarray):
    '''
    Gets the deviation from the annual mean at each step of a time series
    '''
    ### Yearly averages ###
    yearlyMn  = np.zeros(len(timearray))
    yrs       = np.unique(np.floor(timearray))
    nyr       = len(yrs)
    for yy in range(nyr):
        yrinds           = np.where(np.floor(timearray)==yrs[yy])[0]
        yearlyMn[yrinds] = np.mean(valsarray[yrinds])

    ### Calculate monthly effects ###
    fulldevs  = valsarray-yearlyMn
    return(fulldevs)

def fromMonthlyToAnnual(timearray,valsarray):
    '''
    Averages values annualy (only over years with values every month)
    '''
    yrinds0 = [np.where(np.floor(timearray)==tm)[0] for tm in np.unique(np.floor(timearray))]
    yrinds  = [ls for ls in yrinds0 if len(ls)==12]
    if(len(yrinds)==0):
        print('No full year in timearray, need at least one 12-month year')
    timeannual = np.array([0.5+timearray[inds[0]] for inds in yrinds])
    valsannual = np.array([np.mean(valsarray[inds]) for inds in yrinds])
    return(timeannual,valsannual)

def piecewisePolyFitSingleBreak(xvals,yvals,Order,BreakPoint=np.nan):
    '''
    Fit a piecewise polynomial of order Order with one breakpoint
    '''
    MinValPerPeriod = 5 #minimum indices to include in pre and post period to keep breakpoint
    ### Use pwlf to find the breakpoint ###
    pwlf_model   = pwlf.PiecewiseLinFit(xvals,yvals,degree=Order)
    if(np.isnan(BreakPoint)):
        pwlf_model.fit(n_segments=2)
    else:
        pwlf_model.fit_with_breaks([xvals[0],BreakPoint,xvals[-1]])
    training  = pwlf_model.predict(xvals) #predict values from pwlf
    brk       = pwlf_model.fit_breaks[1] #the brakpoint
    iipre     = xvals<=brk
    iipost    = xvals>brk
    if(np.sum(iipre)>=MinValPerPeriod and np.sum(iipost)>=MinValPerPeriod):
        # Enough data points pre- and post-break #
        b0pre     = pwlf_model.predict(xvals[0])  #intercept pre
        b0post    = pwlf_model.predict(brk) #intercept post
        # Fit numpy to training (pwlf beta parameters are not reliable for degree>1) #
        polypre   = np.polyfit(xvals[iipre]-xvals[0],training[iipre]-b0pre,deg=Order)
        polypost  = np.polyfit(xvals[iipost]-brk,training[iipost]-b0post,deg=Order)
    else:
        # Not enough data points pre- and post-break: use a single polynomial #
        brk    = xvals[0]
        iipre  = np.ones(len(xvals)).astype(bool)
        iipost = np.zeros(len(xvals)).astype(bool)
        # Fit with a single segment #
        pwlf_model   = pwlf.PiecewiseLinFit(xvals,yvals,degree=Order)
        pwlf_model.fit(n_segments=1)
        training  = pwlf_model.predict(xvals) #predict values from pwlf
        b0pre     = pwlf_model.predict(xvals[0])  #intercept pre
        polypre   = np.polyfit(xvals-xvals[0],training-b0pre,deg=Order)
        polypost  = np.zeros(Order+1) #no post parameters
        b0post    = 0 #no post parameters

    params    = np.append(np.append(polypre[0:-1],b0pre),
                          np.append(polypost[0:-1],b0post))

    return([params,brk,iipre,iipost])


def evaluateARMAmodel(timearray,valsarray,alphalevel,parray,qarray,trendType='ct',printSummary=False,method='statespace'):
    '''
    Evaluates Bayesian Information Criterion for ARMA models with different
    parameter combinations of AR and MA components
    And calculates if all parameters included in the fit are significant
        (level of significance determined by alphalevel)
    '''

    # Convert single values to 1-element array if needed
    if(len(np.shape(parray))==0):
        parray = np.array([parray])
    if(len(np.shape(qarray))==0):
        qarray = np.array([qarray])
    if(len(parray)!=len(qarray)):
        print('Error in evaluateBICofARMA: parray and qarray should be of same length')
        sys.exit()
    ### Prepare array of BIC values and of number of non-singificant parameters ###
    bicarray    = np.zeros(len(parray)) #prepare array of BIC values
    allsigarray = np.zeros(len(parray)) #prepare binary array evaluating if all params included are significant
    ### Convert to Pandas ###
    vlPandas  = pd.Series(valsarray)
    ### Fit AR models ###
    for ii in range(len(parray)):
        model = ARIMA(vlPandas,order=(parray[ii],0,qarray[ii]),trend=trendType)
        resmodel = model.fit(method=method)
        if(printSummary):
            print(resmodel.summary())
        # Save BIC value and significance of parameters #
        if(np.all(np.array(resmodel.pvalues<=alphalevel))):
            allsigarray[ii] = 1
        bicarray[ii]  = resmodel.bic
    return([bicarray,allsigarray])

def fitAndGenerateARMA(timearray,valsarray,porder,qorder,nsamples,trendType='ct',printSummary=False,method='statespace'):
    '''
    Fits an ARMA model
    AR order is porder
    MA order is qorder
    Generates nsamples random realizations
    Returns a dictionary with:
        ARlagparameters,MAlagparameters,const,trend,noisevariance,randomsamples
    '''
    ### Conversion to Pandas series ###
    vlPandas    = pd.Series(valsarray)
    ### Time scaling (Autoreg assumes timestep=1) ###
    scalingtime = 1/np.mean(np.diff(timearray)) #time scaling (Autoreg assumes timestep=1)
    ### Fit ARMA model ###
    model = ARIMA(vlPandas,order=(porder,0,qorder),trend=trendType)
    resmodel = model.fit(method=method)
    if(printSummary):
        print(resmodel.summary())
    paramsARlag = np.zeros(porder)
    paramsMAlag = np.zeros(qorder)
    for ii in range(porder):
        paramsARlag[ii] = resmodel.params[f'ar.L{ii+1}']
    for ii in range(qorder):
        paramsMAlag[ii] = resmodel.params[f'ma.L{ii+1}']
    dicout = {'ARlagparameters':paramsARlag,'MAlagparameters':paramsMAlag}
    constandtrendterm  = np.zeros(len(timearray))
    if('c' in trendType):
        dicout['const']    = resmodel.params['const']
        constandtrendterm += resmodel.params['const']
    else:
        dicout['const']    = 0.0
    if('t' in trendType):
        dicout['trend']    = resmodel.params['x1']
        constandtrendterm += resmodel.params['x1']*(timearray-timearray[0])*scalingtime
    else:
        dicout['trend']    = 0.0
    dicout['noisevariance'] = copy.deepcopy(resmodel.params['sigma2'])  
    dicout['residuals']     = np.copy(resmodel.resid)
    ### Compute random samples ###
    rdmsamples           = np.zeros((nsamples,len(timearray)))
    rdmsamples[:,0:max(porder,qorder)] = valsarray[0:max(porder,qorder)]
    noiseterms           = np.random.normal(0,np.sqrt(resmodel.params['sigma2']),size=(nsamples,len(timearray)))
    for tt,ttm in enumerate(timearray):
        if(tt>=porder and tt>=qorder):
            rdmsamples[:,tt] += constandtrendterm[tt]
            artermsprev = rdmsamples[:,tt-porder:tt][:,::-1] #reverse to multiply tt-1 by lag1 coef
            artermsprev = artermsprev-constandtrendterm[tt-porder:tt][::-1] #statsmodels ARMA constant and trend terms must be subtracted from AR values
            arlagterms  = paramsARlag*artermsprev
            matermsprev = noiseterms[:,tt-qorder:tt][:,::-1] #reverse to multiply tt-1 by lag1 coef
            malagterms  = paramsMAlag*matermsprev
            if(porder==1):
                arlagterms = arlagterms.flatten()
            else:
                arlagterms = np.sum(arlagterms,axis=1)
            if(qorder==1):
                malagterms = malagterms.flatten()
            else:
                malagterms = np.sum(malagterms,axis=1)
            rdmsamples[:,tt] += (arlagterms+malagterms+noiseterms[:,tt])
    dicout['randomsamples'] = copy.deepcopy(rdmsamples)
    return(dicout)

def generateARMAsamples(timearray,const,trend,arparams,maparams,noisevariance,nsamples,noiseseries=np.nan):
    '''
    Generates nsamples random realizations of an ARMA process
    A pre-computed noise series can be passed as input argument (noiseseries)
    '''
    ### Get orders of AR and MA parts ###
    porder = len(arparams)
    qorder = len(maparams)
    if(porder==0):
        arparams = np.array([0])
        porder   = 1
    if(qorder==0):
        maparams = np.array([0])
        qorder   = 1
    constandtrendterm  = const+trend*(timearray-timearray[0])
    ### Compute random samples ###
    rdmsamples           = np.zeros((nsamples,len(timearray)))
    ii0                  = max(porder,qorder) #first index of computable ARMA
    if(np.all(np.isnan(noiseseries)==True)):
        noiseterms       = np.random.normal(0,np.sqrt(noisevariance),size=(nsamples,len(timearray)))
    else:
        if(np.shape(noiseseries)[0]!=nsamples or np.shape(noiseseries)[1]!=len(timearray)):
            print('Error in generateARMAsamples: input noiseries has wrong dimensions')
            sys.exit()
        else:
            noiseterms = copy.deepcopy(noiseseries)
    rdmsamples[:,0:ii0]  = constandtrendterm[0:ii0]+noiseterms[:,0:ii0]
    for tt,ttm in enumerate(timearray):
        if(tt>=ii0):
            rdmsamples[:,tt] += constandtrendterm[tt]
            artermsprev = rdmsamples[:,tt-porder:tt][:,::-1] #reverse to multiply tt-1 by lag1 coef
            artermsprev = artermsprev-constandtrendterm[tt-porder:tt][::-1] #statsmodels ARMA constant and trend terms must be subtracted from AR values
            arlagterms  = arparams*artermsprev
            matermsprev = noiseterms[:,tt-qorder:tt][:,::-1] #reverse to multiply tt-1 by lag1 coef
            malagterms  = maparams*matermsprev
            if(porder==1):
                arlagterms = arlagterms.flatten()
            else:
                arlagterms = np.sum(arlagterms,axis=1)
            if(qorder==1):
                malagterms = malagterms.flatten()
            else:
                malagterms = np.sum(malagterms,axis=1)
            rdmsamples[:,tt] += (arlagterms+malagterms+noiseterms[:,tt])
    return(rdmsamples)

def mapBias01(tfmodel,tfref):
    '''
    Function mapping bias from [0:+inf] to [0:1]
    '''
    val0   = 1/abs(np.mean(tfmodel)-np.mean(tfref))
    valout = val0/(1+val0)
    return(valout)

def qualMatsOliverHolbrook(timemodel,tfmodel,timereference,tfreference,steppingSsn):
    '''
    Computes Mean, Seasonal and Residual quality matrices values for TF time series
    Follows Eqs.(2,3,4) of Oliver and Holbrook (2014)
    Detrending TF series before any computation
    Computes seasonality effect as the average deviation of steps from the yearly mean (each steppingSsn)
    Note that quality of mean uses modified function for bias
    Note the use of (1+rho)/2 mapping of correlation coefficients to deal with rho<0
    Geometric mean of corr and total abs monthly effects for the seasonality criterion
    Geometric mean of corr and sdevagreement for the residuals criterion
    '''
    ### Find overlapping periods ###
    ttMod   = np.logical_and(timemodel>=timereference[0],timemodel<=timereference[-1])
    ttRef   = np.logical_and(timereference>=timemodel[0],timereference<=timemodel[-1])
    # Adjust first or last index if one time series has one more value #
    if(np.sum(ttRef)==np.sum(ttMod)+1):
        if(abs(timereference[ttRef][0]-timemodel[ttMod][0])>abs(timereference[ttRef][-1]-timemodel[ttMod][-1])):
            ttRef[np.where(ttRef)[0][-1]] = False
        else:
            ttRef[np.where(ttRef)[0][0]] = False
    if(np.sum(ttMod)==np.sum(ttRef)+1):
        if(abs(timereference[ttRef][0]-timemodel[ttMod][0])>abs(timereference[ttRef][-1]-timemodel[ttMod][-1])):
            ttMod[np.where(ttMod)[0][-1]] = False
        else:
            ttMod[np.where(ttMod)[0][0]] = False  
    timeMod = timemodel[ttMod]
    timeRef = timereference[ttRef]
    tfMod0  = tfmodel[ttMod]
    tfRef0  = tfreference[ttRef]
    ### Detrending ###
    mylinregress   = np.polyfit(timeMod,tfMod0,deg=1)
    removed        = mylinregress[0]*timeMod+mylinregress[1]-np.mean(tfMod0)
    tfMod          = tfMod0-removed
    mylinregress   = np.polyfit(timeRef,tfRef0,deg=1)
    removed        = mylinregress[0]*timeRef+mylinregress[1]-np.mean(tfRef0)
    tfRef          = tfRef0-removed
    ### Compute quality of Mean ###
    qualMn  = mapBias01(tfMod,tfRef) 
    ### Compute quality of seasonality ###
    tfnomnMod  = tfMod-np.mean(tfMod)
    tfnomnRef  = tfRef-np.mean(tfRef)
    if(steppingSsn==12):
        fullmonthlyeffectsMod = getMonthlyEffects(timeMod,tfnomnMod)
        fullmonthlyeffectsRef = getMonthlyEffects(timeRef,tfnomnRef)
        sumAbsMonthEffMod     = sum(abs(fullmonthlyeffectsMod[0:12]))
        sumAbsMonthEffRef     = sum(abs(fullmonthlyeffectsRef[0:12]))
    else:
        print('steppingSsn value not implemented yet in qualMatsOliverHolbrook')
        sys.exit()
    corrTermSsn = (1+np.corrcoef(fullmonthlyeffectsMod,fullmonthlyeffectsRef)[0,1])/2
    amplTermSsn = min(sumAbsMonthEffMod/sumAbsMonthEffRef,sumAbsMonthEffRef/sumAbsMonthEffMod)
    qualSsn     = (corrTermSsn*amplTermSsn)**(1/2)
    ### Compute quality of residuals ###
    rsdMod      = tfnomnMod-fullmonthlyeffectsMod
    rsdRef      = tfnomnRef-fullmonthlyeffectsRef
    corrTermRsd = (1+np.corrcoef(rsdMod,rsdRef)[0,1])/2
    sdevTermRsd = min(np.std(rsdMod)/np.std(rsdRef),np.std(rsdRef)/np.std(rsdMod))
    qualRsd     = (corrTermRsd*sdevTermRsd)**(1/2)
    
    return([qualMn,qualSsn,qualRsd])


def streMatsOliverHolbrook(timing,tf0,tf1,steppingSsn):
    '''
    Computes Mean, Seasonal and Residual strength matrices values for TF time series
    Follows Eqs.(19,20,21) of Oliver and Holbrook (2014)
    Detrending TF series before any computation
    Computes seasonality effect as the average deviation of steps from the yearly mean (each steppingSsn)
    Note that strength of mean uses modified function for bias
    Note the use of (1+rho)/2 mapping of correlation coefficients to deal with rho<0
    Note that this assumes same time period for tf0 and tf1 because they should be from the same product
    Geometric mean of corr and total abs monthly effects for the seasonality criterion
    Geometric mean of corr and sdevagreement for the residuals criterion
    '''
    ### Consistency check ###
    if(len(timing)!=len(tf0) or len(timing)!=len(tf1)):
        print('Error in streMatsOliverHolbrook(): timing, tf0 and tf1 should be of same length')
        sys.exit()
    tf0wtr = np.copy(tf0)
    tf1wtr = np.copy(tf1)
    ### Detrending ###
    mylinregress   = np.polyfit(timing,tf0wtr,deg=1)
    removed        = mylinregress[0]*timing+mylinregress[1]-np.mean(tf0wtr)
    tf0            = tf0wtr-removed
    mylinregress   = np.polyfit(timing,tf1wtr,deg=1)
    removed        = mylinregress[0]*timing+mylinregress[1]-np.mean(tf1wtr)
    tf1            = tf1wtr-removed
    ### Compute strength of Mean ###
    streMn = mapBias01(tf0,tf1)
    ### Compute strength of seasonality ###
    tfnomn0  = tf0-np.mean(tf0)
    tfnomn1  = tf1-np.mean(tf1)
    if(steppingSsn==12):
        fullmonthlyeffects0 = getMonthlyEffects(timing,tf0)
        fullmonthlyeffects1 = getMonthlyEffects(timing,tf1)
        sumAbsMonthEff0     = sum(abs(fullmonthlyeffects0[0:12]))
        sumAbsMonthEff1     = sum(abs(fullmonthlyeffects1[0:12]))
    else:
        print('steppingSsn value not implemented yet in streMatsOliverHolbrook')
        sys.exit()
    corrTermSsn = (1+np.corrcoef(fullmonthlyeffects0,fullmonthlyeffects1)[0,1])/2
    amplTermSsn = min(sumAbsMonthEff0/sumAbsMonthEff1,sumAbsMonthEff1/sumAbsMonthEff0)
    streSsn     = (corrTermSsn*amplTermSsn)**(1/2)
    ### Compute strength of residuals ###
    rsd0        = tfnomn0-fullmonthlyeffects0
    rsd1        = tfnomn1-fullmonthlyeffects1
    corrTermRsd = (1+np.corrcoef(rsd0,rsd1)[0,1])/2
    sdevTermRsd = min(np.std(rsd0)/np.std(rsd1),np.std(rsd1)/np.std(rsd0))
    streRsd     = (corrTermRsd*sdevTermRsd)**(1/2)
    
    return([streMn,streSsn,streRsd])


def locMatOliverHolbrook(lat0,lon0,lat1,lon1,DecayScale):
    '''
    Computes Localization matrix values for two pairs of lat,lon coords
    Follows Eq.(22) of Oliver and Holbrook (2014)
    Decay scale should be provided in [m]
    '''
    ### Get coordinates in (x,y) plane ###
    (xx0,yy0) = transformer1.transform(lat0,lon0)
    (xx1,yy1) = transformer1.transform(lat1,lon1)
    ### Compute squared distance ###
    dist2     = (xx0-xx1)**2+(yy0-yy1)**2
    ### Compute localization matrix value ###
    locVal    = np.exp(-1*dist2/(2*DecayScale**2))
    
    return(locVal)

def extrapolateOliverHolbrook(tmodel,tfmeanmodOffsh,tfssnmodOffsh,tfrsdmodOffsh,tref,tfmeanrefOffsh,tfssnrefOffsh,tfrsdrefOffsh,tfrefInsh,steppingSsn,detrendResid,retrieveParams=False):
    '''
    Computes Stat relations between tfrefOffsh and tfrefInsh
    Detrending of time series used for residuals
    Computes the relations for Mean, Seasonality, and Residuals using different gridpoints
    Applies Stat relations to tfmodel
    Using seasonality effect as the average deviation of steps from the yearly mean (each steppingSsn)
    Uses a different extrapolation model for the residuals: calibrating StdDev and values from CDF
        rsd approach: overestimates Offsh-Insh correlation, but: avoids underestim of rsdInsh if low correlation Insh-Offsh in ref
    Boolean option to detrend the residuals prior to deriving rsd relation and applying rsd relation
    Follows Eqs.(7-16) of Oliver and Holbrook (2014)
    '''
    ### Stat relation for Mean ###
    prma = np.mean(tfrefInsh)/np.mean(tfmeanrefOffsh) #see Eq.(7)
    mnInshMod = prma*np.mean(tfmeanmodOffsh) #see Eq.(9)
    if(detrendResid):
        # Detrend TF #
        tfrefInsh0      = np.copy(tfrefInsh)
        mylinregress    = np.polyfit(tref,tfrefInsh0,deg=1)
        removedrefInsh  = mylinregress[0]*tref+mylinregress[1]-np.mean(tfrefInsh0)
        tfrefInsh       = tfrefInsh0-removedrefInsh
    else:
        # Keep tfrefInsh intact #
        pass
    
    ntot  = len(tmodel)
    ### Get seasonal signals from ssn TF time series ###
    if(steppingSsn!=12):
        print('steppingSsn value not implemented yet in streMatsOliverHolbrook')
        sys.exit()
    tfnomnRefInsh    = tfrefInsh-np.mean(tfrefInsh)
    fullssneffModOffsh   = getMonthlyEffects(tmodel,tfssnmodOffsh-np.mean(tfssnmodOffsh)) 
    fullssneffRefOffsh   = getMonthlyEffects(tref,tfssnrefOffsh-np.mean(tfssnrefOffsh)) 
    fullssneffRefInsh    = getMonthlyEffects(tref,tfnomnRefInsh) 
    ### Scaling by sum of abs avoids problem of single effect close to 0 ###
    totssneffRefOffsh    = np.sum(abs(fullssneffRefOffsh[0:12]))
    totssneffRefInsh     = np.sum(abs(fullssneffRefInsh[0:12]))
    prmgamma             = totssneffRefInsh/totssneffRefOffsh
    ### Compute seasonal component of Inshore Model ###
    ssnInshMod           = prmgamma*fullssneffModOffsh
    ### Get residual signals from rsd TF time series ###
    if(detrendResid):
        # Detrend TF #
        mylinregress    = np.polyfit(tmodel,tfrsdmodOffsh,deg=1)
        removedMod      = mylinregress[0]*tmodel+mylinregress[1]-np.mean(tfrsdmodOffsh)
        tf0             = tfrsdmodOffsh-removedMod
    else:
        tf0             = np.copy(tfrsdmodOffsh)
    tf1                 = tf0-np.mean(tf0)
    tfrsdonlymodOffsh0  = tf1-getMonthlyEffects(tmodel,tf1)
    if(detrendResid):
        # Detrend TF #
        mylinregress    = np.polyfit(tref,tfrsdrefOffsh,deg=1)
        removedrefOffsh = mylinregress[0]*tref+mylinregress[1]-np.mean(tfrsdrefOffsh)
        tf0             = tfrsdrefOffsh-removedrefOffsh
    else:
        tf0             = np.copy(tfrsdrefOffsh)
    tf1                 = tf0-np.mean(tf0)
    tfrsdonlyrefOffsh0  = tf1-getMonthlyEffects(tref,tf1)    
    tfrsdonlyrefInsh0   = tfnomnRefInsh-getMonthlyEffects(tref,tfnomnRefInsh)
    sigrsdrefOffsh      = np.std(tfrsdonlyrefOffsh0)
    sigrsdrefInsh       = np.std(tfrsdonlyrefInsh0)
    betaparam           = sigrsdrefInsh/sigrsdrefOffsh
    # Multiplying by the StdDev ratio scales all the values to the chosen normal distri (if mean==0) #
    rsdInshMod0         = betaparam*tfrsdonlymodOffsh0
    if(detrendResid):
        # Add back the trend #
        rsdInshMod          = rsdInshMod0+removedMod
        tfrsdonlymodOffsh   = tfrsdonlymodOffsh0+removedMod
        tfrsdonlyrefOffsh   = tfrsdonlyrefOffsh0+removedrefOffsh
        tfrsdonlyrefInsh    = tfrsdonlyrefInsh0+removedrefInsh

    else:
        # Nothing to be done #
        rsdInshMod          = np.copy(rsdInshMod0) 
        tfrsdonlymodOffsh   = np.copy(tfrsdonlymodOffsh0)
        tfrsdonlyrefOffsh   = np.copy(tfrsdonlyrefOffsh0)
        tfrsdonlyrefInsh    = np.copy(tfrsdonlyrefInsh0)
        
    ### Prepare output ###
    matoutputMod    = np.row_stack((np.mean(tfmeanmodOffsh)*np.ones(ntot),fullssneffModOffsh,tfrsdonlymodOffsh,
                                  mnInshMod*np.ones(ntot),ssnInshMod,rsdInshMod))
    matoutputRef    = np.row_stack((np.mean(tfmeanrefOffsh)*np.ones(len(tfrsdonlyrefOffsh)),fullssneffRefOffsh,tfrsdonlyrefOffsh,
                                    np.mean(tfrefInsh)*np.ones(len(tfrsdonlyrefInsh)),fullssneffRefInsh,tfrsdonlyrefInsh))
    fulloutput    = [matoutputMod,matoutputRef]
    if retrieveParams:
        # Add array of prma,prmgammas,prmdeltas,betaparam to the output list #
        extrapolparams = np.array([prma,prmgamma,betaparam])
        inshparams     = np.concatenate(([mnInshMod],ssnInshMod[0:steppingSsn],
                                          [betaparam*np.std(tfrsdonlymodOffsh)]))
        fulloutput.append(extrapolparams)
        fulloutput.append(inshparams)

    return(fulloutput)



