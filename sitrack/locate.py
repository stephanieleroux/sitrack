### TO DO:
### * use `TheCell` also in `FindContainingCell`
### * their should be one:
###   => FindContainingCell_geo (currently FCC) and FindContainingCell_crt (currently FindContainingCell)

import numpy as np

Pi = 3.141592653589793
EarthRad = 6300. ;
toRad = Pi/180.


def find_ji_of_min(x):
    '''
    # Yes, reinventing the wheel here, but it turns out
    # it is faster this way!
    '''
    k = x.argmin()
    nx = x.shape[1]
    return k//nx, k%nx


def NorthStereoProj( pphi, plam, lam0=0., phi0=90. ):
    #
    # https://mathworld.wolfram.com/StereographicProjection.html
    #
    # to km (since Earth radius provided in km)!!!
    #
    
    zsinPhi0  = np.sin(toRad*phi0)
    zsinPhi   = np.sin(toRad*pphi)
    zcosPhi0  = np.cos(toRad*phi0)
    zcosPhi   = np.cos(toRad*pphi)    
    zcosLambL = np.cos(toRad*(plam - lam0))
    
    zk = 2.*EarthRad / ( 1. + zsinPhi0*zsinPhi + zcosPhi0*zcosPhi*zcosLambL )

    zx = zk * zcosPhi*np.sin(toRad*(plam - lam0))

    zy = zk * ( zcosPhi0*zsinPhi - zsinPhi0*zcosPhi*zcosLambL )


    return zy, zx





def IsInsideQuadrangle( y, x, quad ):
    '''
        Meant to be used with Cartesian coordinates ([m] or [km])

        * quad: shape(4,2), like this: [ [y0,x0], [y1,x1], [y2,x2], [y3,x3] ]

    '''
    n = len(quad)
    if len(quad) != 4:
        print('ERROR: `len(quad) !=: 4`') ; exit(0)
    #
    lInside = False
    p2x = 0.0
    z2y = 0.0
    xints = 0.0
    [z1y,z1x] = quad[0,:]
    #
    for i in range(n+1):
        z2y,z2x = quad[i%n,:]
        if y > min(z1y,z2y):
            if y <= max(z1y,z2y):
                if x <= max(z1x,z2x):
                    if z1y != z2y:
                        xints = (y-z1y)*(z2x-z1x)/(z2y-z1y) + z1x
                    if z1x == z2x or x <= xints:
                        lInside = not lInside
        #
        z1y, z1x = z2y, z2x
        #
    return lInside


def TheCell( pyx, kjiT, pYf, pXf, iverbose=0 ):
    ''' 
        # LOLO: I want it to use Lat,Lon rather than Y,X !!!

        Based on the nearest T-point, find the cell/mesh (defined as the polygon that
        joins 4 neighbor F-points) that contains the target point.
        
        Mind: the indices of the T-point we return is that of the center of the
              identified cell! And in very unusual cases, this is not the same
              as the nearest T-point provided as an input...
       Input:
        * pyx      : (y,x) target point cartesian coordinates [km]
        * kjiT     : (j,i) indices of nearest T-point to (y,x) that was found...
        * pYf, pXf : 2D arrays of cartesian coordinates of the model F-point grid
       Output:
        * lPin    : success or not (boolean)
        * [jT,iT] : actual T-point at the center of the identified cell
        * [[jf1,jf2,jf3,jf4],[if1,if2,if3,if4]]: the 4 F-points that make the vertices
                                                 of the identified cell (ConterClockwize starting from BLC)
    '''
    (zy,zx) = pyx
    (kj,ki) = kjiT
    #
    lPin = False
    kp=0             ; # pass counter
    while (not lPin) and (kp<5):
        kp = kp + 1
        if iverbose>0 and kp>1: print('  * [SeedInit()] search for proper F-cell => test option #'+str(kp)+'!')
        if   kp==1:
            jT,iT =kj,ki   ; # Normal case!!! 99% of all cases !!!
        elif kp==2:
            jT,iT =kj,ki+1 ; # maybe F-cell is the one next to the right?
        elif kp==3:
            jT,iT =kj+1,ki ; # maybe F-cell is the one next above?
        elif kp==4:
            jT,iT =kj,ki-1 ; # maybe F-cell is the one next to the left?
        elif kp==5:
            jT,iT =kj-1,ki ; # maybe F-cell is the one next below?
        #
        # Based on this @center T-point (here `jT,iT`), the proper cell/mesh
        # should be the polygon defined by the 4 following F-points:
        # (indexing is anti-clockwize, starting from bottom left F-point)
        [jf1,jf2,jf3,jf4] = [ jT-1, jT-1, jT, jT  ]
        [if1,if2,if3,if4] = [ iT-1, iT,   iT, iT-1 ]                    
        #
        zquad = np.array( [ [pYf[jf1,if1],pXf[jf1,if1]], [pYf[jf2,if2],pXf[jf2,if2]],
                            [pYf[jf3,if3],pXf[jf3,if3]], [pYf[jf4,if4],pXf[jf4,if4]] ] )
        lPin = IsInsideQuadrangle( zy, zx, zquad )
        #
    if iverbose>0:
        if kp>1 and lPin: print('        => option #'+str(kp)+' did work!  :)')
    #
    return lPin, [jT,iT], np.array( [ [jf1,if1],[jf2,if2],[jf3,if3],[jf4,if4] ], dtype=int )


# This is to become a generic "Find containing Cell" for NEMO...
# => can be used to find the T-centric cell (with 4 vertices being F-points)
# => or to find the F-centric cell (with 4 vertices being T-points) => what we need for linear interpolation of a T-field onto point inside this cell...
def FCC( pntGcoor, pLat, pLon, pLatC, pLonC, cellType='T', rd_found_km=10., resolkm=[],
         ji_prv=(), np_box_r=10, max_itr=5, pntID=None, iverbose=0 ):
    '''
        Provided an input target-point geographic coordinates, locate the grid mesh cell that
        includes this point.
        The grid mesh is of type Arakawa C-grid, typically a NEMO/ORCA type of grid...

    INPUT:
             * pntGcoor  : GPS coordinates (lat,lon) of target point    ([real],[real])
             * pLat      : array of source grid lat.  @T/F-points if `cellType='T/F'`) [real]
             * pLon      : array of source grid long. @T/F-points if `cellType='T/F'`) [real]
             * pLatC     : array of source grid lat.  @F/T-points if `cellType='T/F'`) [real]
             * pLonC     : array of source grid long. @F/T-points if `cellType='T/F'`) [real]
             * resolkm   : array of source grid approximate local resolution [km] [real]
                           because grids like ORCA of NEMO can have strong spatial heterogenity of resolution...
    RETURNS:
             * j,i : indices of the grid mesh cell containing the target point
                     => -1,-1 if something went wrong or point was not found

    TODO: must take care of boundaries when doing local square box!

    '''
    nhb = 2 ; # half width of local box in number of points...
    #
    (zlat,zlon) = pntGcoor
    
    lbla = ( iverbose>0 and pntID )
    
    cP0 = cellType      ; # string for type of center point 
    if cellType=='T':
        cP4C = 'F'      ; # string for type of corner/vertices points
    elif cellType=='F':
        cP4C = 'T'      ; # string for type of corner/vertcices points
    else:
        print('ERROR [FCC]: for now we just expect the mesh to be centered on "T" or "F" points.')
        exit(0)

    icancel = 0
    
    # First, find the nearest `cellType`-point to `pntGcoor`:
    (jX, iX) = NearestPoint( pntGcoor, pLat, pLon, rd_found_km=rd_found_km, resolkm=resolkm,
                             ji_prv=ji_prv, np_box_r=np_box_r, max_itr=max_itr )

    if jX>=0 and iX>=0:
        # Ok a nearest point was found!                
        if iverbose>0:
            print('     ==> nearest '+cP0+'-point for ',zlat,zlon,' on target grid:', jX, iX, '==> lat,lon:',
                  round(pLat[jX,iX],3), round(pLon[jX,iX],3))


        # To avoid having problem working with `lat,lon` we apply a polar stereographic projection to
        # a small square domain that surround the nearest point found:
        #  * lat,lon => y,k (km)        
        zYC, zXC = NorthStereoProj( pLatC[jX-nhb:jX+3,iX-nhb:iX+3], pLonC[jX-nhb:jX+3,iX-nhb:iX+3],
                                        lam0=pLon[jX,iX], phi0=pLat[jX,iX] )

        # Projected target point:
        zy0, zx0 = NorthStereoProj( np.array([zlat]) , np.array([zlon]),
                                    lam0=pLon[jX,iX], phi0=pLat[jX,iX] )
            
        lPin, [jM,iM], KVRTCS = TheCell( (zy0,zx0), (nhb,nhb), zYC, zXC, iverbose=2 )
        
        # Back to the global frame
        jX = jM + jX - nhb
        iX = iM + iX - nhb
        KVRTCS[:,0] = KVRTCS[:,0] + jX - nhb
        KVRTCS[:,1] = KVRTCS[:,1] + iX - nhb
        
        if not lPin:
            print('WARNING [SeedInit()]: could not find the proper F-point cell!!!')
            print('         => when lookin for point:',zlat,zlon)
            icancel = 0

    else:
        if iverbose>0:
            print('  [FCC] ==> NO nearest '+cP0+'-point found for ',zlat,zlon,'!')


    #
    return [jX,iX], KVRTCS
    


def NearestPoint( pntGcoor, pLat, pLon, rd_found_km=10., resolkm=[], ji_prv=(), np_box_r=10, max_itr=5 ):
    '''
    # * pntGcoor : GPS coordinates (lat,lon) of target point    ([real],[real])
    # * pLat        : array of source grid latitude            2D numpy.array [real]
    # * pLon        : array of source grid longitude           2D numpy.array [real]
    # * resolkm   : array of source grid approximate local resolution [km] 2D numpy.array [real]
    #               because grids like ORCA of NEMO can have strong spatial heterogenity of resolution...
    '''
    from .util import Haversine
    #
    (Ny,Nx) = pLat.shape
    if np.shape(pLon) != (Ny,Nx):
        print('ERROR [NearestPoint]: `pLat` & `pLon` do not have the same shape!')
        exit(0)
    l2Dresol = ( np.shape(resolkm)==(Ny,Nx) )
    lbox     = ( len(ji_prv)==2 )
    #
    (latP,lonP) = pntGcoor
    #        
    if lbox:
        (j_prv,i_prv) = ji_prv
        j1, j2 = max(j_prv-np_box_r,0), min(j_prv+np_box_r+1,Ny)
        i1, i2 = max(i_prv-np_box_r,0), min(i_prv+np_box_r+1,Nx)
    else:
        (j1,i1 , j2,i2) = (0,0 , Ny,Nx)
    #
    jy, jx = -1,-1 ; # "not found" flag value...
    lfound = False    
    rfnd   = rd_found_km
    igo    = 0
    #
    while (not lfound) and igo<max_itr :
        igo = igo + 1
        if lbox and igo>1:
            (j1,i1 , j2,i2) = (0,0 , Ny,Nx) ; # Falling back on whole domain for second pass...
        xd = Haversine( latP, lonP,  pLat[j1:j2,i1:i2], pLon[j1:j2,i1:i2] )
        jy, jx = find_ji_of_min( xd )
        #
        if igo==1 and l2Dresol: rfnd = 0.5*resolkm[jy,jx]
        #
        if not lbox and igo==1: igo=2 ; # we jump one round because no need of the pass to global...
        #
        lfound = ( xd[jy,jx] < rfnd )
        if igo>1 and not lfound:
            rfnd = 1.2*rfnd ; # increasing validation distance criterion by 20 %
    #
    if lbox:
        jy, jx = jy+j1, jx+i1 ; # found in the zoom box => translate to indices in whole domain:
    #
    if jy<0 or jx<0 or jy>=Ny or jx>=Nx or igo==max_itr:
        print('    WARNING [NearestPoint()]: did not find a nearest point for target point ',latP,lonP,' !')
        print('            => last tested distance criterions =', rfnd,' km')
        jy, jx = -1,-1
    #
    return (jy, jx)



def FindContainingCell( pyx, kjiT, pYf, pXf, iverbose=0 ):
    ''' 
        Based on the nearest T-point, find the cell/mesh (defined as the polygon that
        joins 4 neighbor F-points) that contains the target point.
        
        Mind: the indices of the T-point we return is that of the center of the
              identified cell! And in very unusual cases, this is not the same
              as the nearest T-point provided as an input...
       Input:
        * pyx      : (y,x) target point cartesian coordinates [km]
        * kjiT     : (j,i) indices of nearest T-point to (y,x) that was found...
        * pYf, pXf : 2D arrays of cartesian coordinates of the model F-point grid
       Output:
        * lPin    : success or not (boolean)
        * [jT,iT] : actual T-point at the center of the identified cell
        * [[jf1,jf2,jf3,jf4],[if1,if2,if3,if4]]: the 4 F-points that make the vertices
                                                 of the identified cell (ConterClockwize starting from BLC)
    '''
    (zy,zx) = pyx
    (kj,ki) = kjiT
    #
    lPin = False
    kp=0             ; # pass counter
    while (not lPin) and (kp<5):
        kp = kp + 1
        if iverbose>0 and kp>1: print('  * [SeedInit()] search for proper F-cell => test option #'+str(kp)+'!')
        if   kp==1:
            jT,iT =kj,ki   ; # Normal case!!! 99% of all cases !!!
        elif kp==2:
            jT,iT =kj,ki+1 ; # maybe F-cell is the one next to the right?
        elif kp==3:
            jT,iT =kj+1,ki ; # maybe F-cell is the one next above?
        elif kp==4:
            jT,iT =kj,ki-1 ; # maybe F-cell is the one next to the left?
        elif kp==5:
            jT,iT =kj-1,ki ; # maybe F-cell is the one next below?
        #
        # Based on this @center T-point (here `jT,iT`), the proper cell/mesh
        # should be the polygon defined by the 4 following F-points:
        # (indexing is anti-clockwize, starting from bottom left F-point)
        [jf1,jf2,jf3,jf4] = [ jT-1, jT-1, jT, jT  ]
        [if1,if2,if3,if4] = [ iT-1, iT,   iT, iT-1 ]                    
        #
        zquad = np.array( [ [pYf[jf1,if1],pXf[jf1,if1]], [pYf[jf2,if2],pXf[jf2,if2]],
                            [pYf[jf3,if3],pXf[jf3,if3]], [pYf[jf4,if4],pXf[jf4,if4]] ] )
        lPin = IsInsideQuadrangle( zy, zx, zquad )
        #
    if iverbose>0:
        if kp>1 and lPin: print('        => option #'+str(kp)+' did work!  :)')
    #
    return lPin, [jT,iT], [[jf1,jf2,jf3,jf4],[if1,if2,if3,if4]]

