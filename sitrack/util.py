#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import exit
#from os import path ; #, mkdir
import numpy as np


def chck4f( cf ):
    from os.path import exists
    if not exists(cf):
        print(' ERROR [chck4f()]: file '+cf+' does not exist!')
        exit(0)

#K
def epoch2clock( it, precision='s' ):
    from datetime import datetime as dt
    from datetime import timezone
    #
    it = int(it)
    if   precision=='s':
        ct = dt.fromtimestamp(it, timezone.utc).strftime("%Y-%m-%d_%H:%M:%S")
    elif precision=='m':
        ct = dt.fromtimestamp(it, timezone.utc).strftime("%Y-%m-%d_%H:%M")
    elif precision=='h':
        ct = dt.fromtimestamp(it, timezone.utc).strftime("%Y-%m-%d_%H")
    elif precision=='D':
        ct = dt.fromtimestamp(it, timezone.utc).strftime("%Y-%m-%d")
    else:
        print('ERROR [epoch2clock]: unknown precision "'+precision+'" !')
        exit(0)
    return str(ct)
#K
def clock2epoch( cdate ):
    from datetime import datetime as dt
    from datetime import timezone
    it = dt.strptime(cdate, "%Y-%m-%d_%H:%M:%S").replace(tzinfo=timezone.utc)
    return int(it.timestamp())

        
def degE_to_degWE( X ):
    '''
    # From longitude in 0 -- 360 frame to -180 -- +180 frame...
    '''
    if np.shape( X ) == ():
        # X is a scalar
        from math import copysign
        return     copysign(1., 180.-X)*        min(X,     abs(X-360.))
    else:
        # X is an array
        return np.copysign(1., 180.-X)*np.minimum(X, np.abs(X-360.))


def DateString( date1, date2, dformat='YYYY-MM-DD_hh:mm:00', returnShort=False ):
    '''
        * date1, date2 : 2 dates in either of these 2 formats `YYYYMMDD` or `YYYYMMDD_hh:mm` !
    '''
    ldp1, lhp1 = len(date1)==8, len(date1)==14
    ldp2, lhp2 = len(date2)==8, len(date2)==14

    cY1,  cY2  = date1[0:4], date2[0:4]
    cmm1, cmm2 = date1[4:6], date2[4:6]
    cdd1, cdd2 = date1[6:8], date2[6:8]

    chhmm1, chhmm2 = '00:00', '00:00'
    if lhp1: chhmm1 = date1[9:11]+':'+date1[12:14]
    if lhp2: chhmm2 = date2[9:11]+':'+date2[12:14]

    cdt1 = str.replace(dformat,'YYYY-MM-DD',cY1+'-'+cmm1+'-'+cdd1)
    cdt2 = str.replace(dformat,'YYYY-MM-DD',cY2+'-'+cmm2+'-'+cdd2)
    cdtL1 = str.replace(cdt1,'hh:mm',chhmm1) ; # ex: 2001-04-22_02:30
    cdtL2 = str.replace(cdt2,'hh:mm',chhmm2)

    cdtS1 = cY1+cmm1+cdd1 ; # ex: 20010422
    cdtS2 = cY2+cmm2+cdd2

    if returnShort:
        return cdtL1, cdtL2,  cdtS1, cdtS2
    else:
        return cdtL1, cdtL2
    
    
def Haversine( plat, plon, xlat, xlon ):
    '''
    # ! VECTOR VERSION !
    # Returns the distance in km at the surface of the earth
    # between two GPS points (degreesN, degreesE)
    # (plat,plon)  : a point
    # xlat, xlon : 2D arrays
    #
    # Here we do not need accuracy on Earth radius, since the call
    # to this function is suposely made to find nearest point
    '''
    to_rad = 3.141592653589793/180.
    R = 6360. ; # radius of Earth in km (sufficient to take it as a constant here...)
    #
    a1 = np.sin( 0.5 * ((xlat - plat)*to_rad) )
    a2 = np.sin( 0.5 * ((xlon - plon)*to_rad) )
    a3 = np.cos( xlat*to_rad ) * np.cos(plat*to_rad)
    #
    return 2.*R*np.arcsin( np.sqrt( a1*a1 + a3 * a2*a2 ) )



def Dist2Coast( rlon, rlat, plon, plat, pdist2coat ):
    '''
       Returns the distance to the nearest coast of a given point (lon,lat)
        INPUT:
          * rlon, rlat: coordinates (scalars) [degrees East], [degrees West]
          * plon:       1D array of longitudes assosiated to `pdist2coat`
          * plat:       1D array of latitudes  assosiated to `pdist2coat`
          * pdist2coat: 2D array containing rasterized distance to coast (originally read in NC file) [km]
    '''
    rx = np.mod( rlon, 360. ) ; # [0:360] frame
    vx = np.mod( plon, 360. ) ; # [0:360] frame
    # Are we dangerously close to the [0--360] cut?
    #  => then will work in the [-180:180] frame:
    if rx > 355.:
        rx = degE_to_degWE( rx )
        vx = degE_to_degWE( vx )
        #print(' rlon, rlat =', rx, rlat)
    ip = np.argmin( np.abs(  vx[:] - rx  ) )
    jp = np.argmin( np.abs(plat[:] - rlat) )
    #print(' ip, jp =', ip, jp)
    #print(' Nearest lon, lat =', plon[ip], plat[jp])
    del vx, rx
    return max( pdist2coat[jp,ip] , 0. )


def TimeBins4Scanning( pdt1, pdt2, pdt,  iverbose=0 ):
    '''
       Need a time axis with bins to look for buoys within...

       Input:
                * pdt1, pdt2 : start & end time ([s] UNIX epoch time)
                * pdt        : width of bin ([s] UNIX epoch time)
       Returns:
                * nB :  number of bins
                * vTB : array(nB,3) axe==0 => time at center of bin      ([s] UNIX epoch time)
                *                        axe==1 => time at lower bound of bin ([s] UNIX epoch time)
                *                        axe==2 => time at upper bound of bin ([s] UNIX epoch time)
    '''
    #
    zhdt = pdt/2.
    #
    nB = int(round((pdt2 - pdt1) / pdt))
    print('\n *** New fixed time axis to use to scan data:\n    ===> nB = '+str(nB)+' time bins!')
    #
    vTB = np.zeros((nB,3), dtype=int  ) ; # `*,0` => precise time | `*,1` => bound below | `*,2` => bound above
    #
    vTB[0,0] =  pdt1 + zhdt    ; # time at center of time bin
    for jt in range(1,nB):
        tt = vTB[jt-1,0] + pdt
        vTB[jt,0] = tt                ; # time at center of time bin
    # Time bins bounds:
    vTB[:,1] = vTB[:,0] - zhdt
    vTB[:,2] = vTB[:,0] + zhdt
    #
    if iverbose>0:
        for jt in range(nB):
            print(" --- jt="+'%3.3i'%(jt)+": * center of bin => ",vTB[jt,0]," => ",epoch2clock(vTB[jt,0]))
            print("             * bin bounds    => "+epoch2clock(vTB[jt,1])+" - "+epoch2clock(vTB[jt,2])+"\n")
    #
    return nB, vTB



def OrderCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''

    # sort the points based on their x-coordinates
    isortX  = np.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]

    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]

    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = np.argsort(leftMost[:,1])
    leftMost   = leftMost[isortL,:]
    (tl, bl)   = leftMost
    [i1l, i2l] = isortML[isortL]

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = np.argsort(rghtMost[:,1])
    rghtMost   = rghtMost[isortR,:]
    (tr, br)   = rghtMost
    [i1r, i2r] = isortMR[isortR]

    #idx = np.concatenate([idxL,idxR])
    #isort = np.array([isortX[i] for i in np.concatenate([isortL,isortR])])

    del isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR

    # return the coordinates in top-left, top-right,
    # bottom-right, and bottom-left order
    return np.array([tl, tr, br, bl], dtype="float32") , np.array([i1l, i1r, i2r, i2l])


def OrderCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''
    zz = xcoor.copy()
    iz = np.array([0,1,2,3])
    zt, it = OrderCW(xcoor)
    #
    zz[1:4,:] = zt[:0:-1,:]
    iz[1:4]   = it[:0:-1]
    del zt, it
    return zz, iz



def SortIndicesCCW(xcoor):
    '''
    ## Sort points defined by their [x,y] coordinates counter-clockwise
    ##
    ## https://gist.github.com/flashlib/e8261539915426866ae910d55a3f9959
    ##
    ## INPUT:
    ##        xcoor: 2D array of coordinates of shape (N,2) => N `[x,y]` coordinates
    ##
    '''
    # sort the points based on their x-coordinates
    isortX  = np.argsort(xcoor[:,0])
    xSorted = xcoor[isortX,:]

    # grab the left-most and right-most points from the sorted
    # x-roodinate points
    leftMost = xSorted[:2,:]
    isortML  =  isortX[:2]
    rghtMost = xSorted[2:,:]
    isortMR  =  isortX[2:]

    # now, sort the left-most coordinates according to their
    # y-coordinates so we can grab the top-left and bottom-left
    # points, respectively
    isortL     = np.argsort(leftMost[:,1])
    [i1l, i2l] = isortML[isortL]

    # if use Euclidean distance, it will run in error when the object
    # is trapezoid. So we should use the same simple y-coordinates order method.

    # now, sort the right-most coordinates according to their
    # y-coordinates so we can grab the top-right and bottom-right
    # points, respectively
    isortR     = np.argsort(rghtMost[:,1])
    [i1r, i2r] = isortMR[isortR]

    #idx = np.concatenate([idxL,idxR])
    #isort = np.array([isortX[i] for i in np.concatenate([isortL,isortR])])

    del xSorted, isortX, isortL, isortR, leftMost, rghtMost, isortML, isortMR

    return np.array([i1l, i1r, i2r, i2l])



def IndXYClones( pXY, rmask_val=-999. ):
    '''
         Input: `pXY` is a 2D array of float containing x,y coordinates for Np points => shape=(Np,2)

             ==> will return the 1D integer vector containing row (points) indices (along Np)
                 corresponding to the location of points that already exist once, so the can be
                 suppresses or masked later on... (only a single occurence of a point coord. is
                 kept). While ignoring masked values flagged with `rmask_val`
    '''
    (_,nc) = np.shape(pXY)
    if (nc !=2 ):
        print('ERROR [IndXYClones]: second dimmension of pXY must be 2! (coordinates)')
        exit(0)

    # Row indices that exclude masked points:
    (IXvalid,_) = np.where( pXY != [rmask_val,rmask_val] )
    IXvalid = IXvalid[::2]

    # Row indices that select masked points:
    (IXmaskd,_) = np.where( pXY == [rmask_val,rmask_val] )
    IXmaskd = IXmaskd[::2]

    _,Iunq = np.unique( pXY, axis=0, return_index=True )

    # keep the values of `Iunq` that are not in `IXmaskd`:
    IXokunq = np.setdiff1d( Iunq, IXmaskd ) ;

    # keep the values of `IXvalid` that are not in `IXokunq`:
    # => that's the indices of points that could be removed to shrink `pXY` !
    return np.setdiff1d( IXvalid, IXokunq )


def GetRidOfXYClones( pXY, rmask_val=-999. ):
    '''
         Input: `pXY` is a 2D array of float containing x,y coordinates for Np points => shape=(Np,2)

             ==> will return updated version of coordinates array `pXY`, rid of points already existing (clones)
    '''
    (_,nc) = np.shape(pXY)
    if (nc !=2 ):
        print('ERROR [GetRidOfXYClones]: second dimmension of pXY must be 2! (coordinates)')
        exit(0)

    # Row indices that select masked points:
    (IXmaskd,_) = np.where( pXY == [rmask_val,rmask_val] )
    IXmaskd = IXmaskd[::2]

    _,Iunq = np.unique( pXY, axis=0, return_index=True )

    # keep the values of `Iunq` that are not in `IXmaskd`:
    IXokunq = np.setdiff1d( Iunq, IXmaskd ) ;

    # keep the values of `IXvalid` that are not in `IXokunq`:
    # => that's the indices of points that could be removed to shrink `pXY` !
    return pXY()



def SubSampCloud( rd_km, pCoor ):
    '''
       * pCoor: [X,Y] !!!
    '''
    from gudhi import subsampling as sbspl
    #
    cerr = 'ERROR [util.SubSampCloud()]: '
    # Sanity check of input:
    if rd_km <= 0. or rd_km > 2000:
        print(cerr+'silly value for `rd_km`:',rd_km)
        exit(0)
    (Nb0,n2) = np.shape(pCoor)
    if (n2 != 2 ):
        print(cerr+'second dimmension of `pCoor` must be 2 !')
        exit(0)
    
    zCoor = np.array( sbspl.sparsify_point_set( pCoor, min_squared_dist=rd_km*rd_km ) )
    (Nb,_) = np.shape(zCoor)

    # Retrieve corresponding indices for selected points:
    idxleft = np.zeros(Nb, dtype=int)
    for i in range(Nb):
        (idx,_) = np.where( pCoor[:,:]==zCoor[i,:] )
        idxleft[i] = idx[0]

    return Nb, zCoor, idxleft
    #if l_do_cgeo:
    #    zgeo = np.array( [ pLonLat[ileft,0], pLonLat[ileft,1] ] ).T
    #    if l_do_name:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft], zgeo, pNames[ileft]
    #    else:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft], zgeo
    #else:
    #    if l_do_name:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft], pNames[ileft]
    #    else:
    #        return Nb, zCoor, pIDs[ileft], ptime[ileft]




def StdDev( pmean, pX ):
    zz = pX[:] - pmean
    return np.sqrt( np.mean( zz*zz ) )





def Geo2CartNPSkm1D( pcoorG, lat0=70., lon0=-45. ):
    '''
         => from Geo coor. (lon,lat)[degrees] to cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj!
    '''
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    (_,n2) = np.shape(pcoorG)
    if n2!=2:
        print(' ERROR [Geo2CartNPSkm1D()]: input array `pcoorG` has a wrong a shape!')
        exit(0)
    #
    print('LOLO Geo2CartNPSkm1D => lat0, lon0 =', lat0, lon0)
    
    crs_src = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    crs_trg = NorthPolarStereo(central_longitude=lon0, true_scale_latitude=lat0) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    #
    zx,zy,_ = crs_trg.transform_points(crs_src, pcoorG[:,1], pcoorG[:,0]).T
    #
    return np.array([ zy/1000., zx/1000. ]).T

#K
def CartNPSkm2Geo1D( pcoorC, lat0=70., lon0=-45. ):
    '''
         => from cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj to Geo coor. (lon,lat)[degrees] !
    '''
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    (_,n2) = np.shape(pcoorC)
    if n2!=2:
        print(' ERROR [CartNPSkm2Geo1D()]: input array `pcoorC` has a wrong a shape!')
        exit(0)    
    #
    crs_src = NorthPolarStereo(central_longitude=lon0, true_scale_latitude=lat0) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    crs_trg = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    #
    zlon,zlat,_ = crs_trg.transform_points(crs_src, 1000.*pcoorC[:,1], 1000.*pcoorC[:,0]).T
    #
    return np.array([zlat, zlon]).T




def ConvertGeo2CartesianNPSkm( plat, plon, lat0=70., lon0=-45. ):
    '''
         => from Geo coor. (lon,lat)[degrees] to cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj!
    '''
    #
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    crs_src = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    crs_trg = NorthPolarStereo(central_longitude=lon0, true_scale_latitude=lat0) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    #
    ndim = len(np.shape(plon))
    #
    zx,zy,_ = crs_trg.transform_points(crs_src, plon, plat).T
    #
    if ndim==2:
        zx,zy = zx.T, zy.T
        #
    return zy/1000., zx/1000.



def ConvertCartesianNPSkm2Geo( pY, pX, lat0=70., lon0=-45. ):
    '''
         => from cartesian (x,y)[km] with RGPS' `NorthPolarStereo` proj to Geo coor. (lon,lat)[degrees] !
    '''
    #
    from cartopy.crs import PlateCarree, NorthPolarStereo
    #
    crs_src = NorthPolarStereo(central_longitude=lon0, true_scale_latitude=lat0) ; # that's (lon,lat) to (x,y) RGPS ! (info from Anton)
    crs_trg = PlateCarree() ;                                                   # this geographic coordinates (lat,lon)
    #
    ndim = len(np.shape(pX))
    #
    zlon,zlat,_ = crs_trg.transform_points(crs_src, 1000.*pX, 1000.*pY).T
    #
    if ndim==2:
        zlon,zlat = zlon.T, zlat.T
        #
    return zlat, zlon



def CheckTimeConsistencyQuads( kF, QD, time_dev_from_mean_allowed, iverbose=0 ):
    '''
        We look at the position DATE of all the points composing a group of quads
        (Quad class `QD`) and return a message error if some points are too far 
        in time from the mean of all the points...

        * kF: file number
        * QD: Quad class loaded from file `kF`
        * time_dev_from_mean_allowed: maximum time deviation from the mean allowed
    '''
    if iverbose>0:
        print('\n *** In file #'+str(kF)+':')
    #
    if len(np.shape(QD.PointTime))>1:
        print('ERROR [CheckTimeConsistencyQuads()]: wrong shape for time array of Quad! (should be 1D!) => ',np.shape(QD.PointTime))
        exit(0)
    if iverbose>0: print('     (time_dev_from_mean_allowed =', time_dev_from_mean_allowed/60.,' minutes)' )
    #for zt in QD.PointTime:  print(epoch2clock(zt)); #lolo
    rTmean = np.mean(QD.PointTime)
    if iverbose>0: cTmean = epoch2clock(rTmean)
    rStdDv = StdDev(rTmean, QD.PointTime)
    if iverbose>0: print('     => Actual Mean time =',cTmean,', Standard Deviation =',rStdDv/60.,' minutes!')
    zadiff = np.abs(QD.PointTime-rTmean)
    zdt = np.max(zadiff)/60.
    if iverbose>0: print('     ==> furthest point is '+str(round(zdt,2))+' minutes away from mean!')
    #
    if np.any(zadiff>time_dev_from_mean_allowed):
        print('ERROR [CheckTimeConsistencyQuads()]: some points read in file #'+str(kF)+' are too far (in time) from mean time!')
        print('                                     => max allowed deviation =',time_dev_from_mean_allowed/3600,'hours')
        (idxFU,) = np.where(zadiff>time_dev_from_mean_allowed)
        dmaxFU = np.max(zadiff[idxFU])
        print('                                     => max found here =',dmaxFU/3600,'hours')
        exit(0)
    else:
        if iverbose>0: print('     => ok! No points further than '+str(time_dev_from_mean_allowed/60.)+' minutes from mean...')
    #
    return rTmean, rStdDv







def CancelTooClose( krec, rdkm, plat, plon, pmsk, NbPass=2, iverbose=0 ):
    '''
       Eliminates points too close to each other (one of the 2, keep the one with longest series!)
       * krec: record to focus on
       * rdkm: distance criterion for elimination [km]
       * plat, plon: 2D arrays of lat,lon [Nrec,Nbuoy]
       * pmsk      : 2D arrays of mask [Nrec,Nbuoy]
    '''
    #
    nB = np.sum(pmsk[krec,:]) ; # Number of valid buoys at specified records
    #
    for jp in range(NbPass):
        print('\n *** Applying initial overlap cleaning at the scale of '+str(rdkm)+' km')
        zbmask = np.zeros(nB, dtype='i1') + 1
        zlat,zlon = plat[krec,:],plon[krec,:]
        zmsk = pmsk.copy()
        #
        for jb in range(nB):
            if jb%1000==0: print('       pass #'+str(jp+1)+'.... '+str(jb)+' / '+str(nB)+'...')
            # All the following work should be done if present buoy has not been cancelled yet
            if zbmask[jb] == 1:
                vdist = Haversine( zlat[jb], zlon[jb], zlat, zlon ) ; # estimates of distances (km) with all other buoy at this same time record:
                #                        
                if vdist[jb] != 0.: print(' PROBLEM: distance with yourself should be 0!'); exit(0)
                vdist[jb] = 9999.     ; # need to mask itself (distance = 0!)
                rdmin = np.min(vdist)
                if iverbose>2: print('    ==> closest neighbor buoy is at '+str(round(rdmin,2))+' km !')
                #
                if rdmin < rdkm:
                    # There is at least 1 buoy too close!
                    # We only deal wit 2 buoys at the time the one we are dealing with and the closest one (otherwize could remove to many)
                    (idx_2c,) = np.where( vdist == rdmin )
                    j2c = idx_2c[0]
                    # Now we have to look at the one of the two that has the longest record:
                    nr1, nr2 = np.sum(zmsk[:,jb]), np.sum(zmsk[:,j2c])
                    if iverbose>2: print(' Nb. of valid records for this buoy and the one too close:',nr1,nr2)
                    if nr1<nr2: j2c=jb ; # we should cancel this one not the found one!
                    zbmask[j2c] = 0
                    zmsk[krec:,j2c] = 0
    ### for jp in range(NbPass)

    nBn = np.sum(zmsk[krec,:])
    print('      => we remove '+str(nB-nBn)+' buoys at all records!')
    (idx_keep,) = np.where(zmsk[krec,:]==1)
    #
    return nBn, idx_keep


