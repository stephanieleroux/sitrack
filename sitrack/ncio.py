
from sys import argv, exit
from os import path ; #, environ, mkdir
import numpy as np

#from re import split
from netCDF4 import Dataset

from .util import chck4f, ConvertGeo2CartesianNPSkm
from .util import epoch2clock as e2c



tunits_default = 'seconds since 1970-01-01 00:00:00'

list_required_var_RGPS = [ 'index', 'x', 'y', 'lon', 'lat', 'time', 'idstream', 'streams' ]

FillValue = -9999.

# Naming convention for Cartesian coordinates in km:
nm_y_t, nm_x_t = 'y_pos_t', 'x_pos_t' ; # at center of grid cell (T-point)
nm_y_f, nm_x_f = 'y_pos_f', 'x_pos_f' ; # at upper-right corner of grid cell (F-point)
nm_y_u, nm_x_u = 'y_pos_u', 'x_pos_u' ; # at center of rhs cell edge (U-point)
nm_y_v, nm_x_v = 'y_pos_v', 'x_pos_v' ; # at center of upper cell edge (V-point)



def GetNameTimeDim( fnc_hndl ):
    '''
    * fnc_hndl: netCDF file handle !
    '''
    list_dim = list(fnc_hndl.dimensions.keys())
    for cdim in [ 'time_counter', 'time', None ]:
        if cdim in list_dim: break
    if cdim: print('    * [GetNameTimeDim]:    => time is called: "'+cdim+'"')
    return cdim


def ConvertTimeToEpoch( vtime, units, calendar ):
    '''
    Convert a time 1D array with unexpected units, such as "days since ..." to a 1D array of UNIX aka Epoch time (seconds since 1970-01-01_00:00:00)
    '''
    from netCDF4  import num2date
    from calendar import timegm
    #
    time_convert = num2date( vtime[:], units, calendar) ; # converts to to python datetime object!
    zepoch = []    
    for tt in time_convert:
        zepoch.append(timegm(tt.timetuple()))
    #
    return np.array( zepoch , dtype='i4' )

    

def GetModelGrid( fNCmeshmask, alsoF=False ):

    l_CartesianCoordPresent = False ; # by default we assume that Cartesian X,Y coordinates are not present in `fNCmeshmask` !
    
    chck4f( fNCmeshmask)
    
    # Reading mesh metrics into mesh-mask file:
    with Dataset(fNCmeshmask) as id_mm:
        #
        # Checking whether Cartesian coordinates are already present or not in `fNCmeshmask`:
        listv = list(id_mm.variables.keys())
        #
        l_CartesianCoordPresent = ( (nm_y_t in listv) and (nm_x_t in listv) and (nm_y_f in listv) and (nm_x_f in listv ) )
        if l_CartesianCoordPresent:
            print('    * [GetModelGrid]: Cartesian coordinates are already present in file "'+fNCmeshmask+'" !')
            print('                      => will skip applying a projection!')
        else:
            print('    * [GetModelGrid]: WARNING! Cartesian coordinates (x,y [km]), not found in "'+fNCmeshmask+'" !')
            print('                      => will apply a projection to create them!')
            print('                      => if it is a mistake and they do exist, then they should be named:')
            print('                         "'+nm_y_t+'", "'+nm_x_t+'", "'+nm_y_f+'", and "'+nm_x_f+'"')
            print('                         and be in km, not m !')
        #
        ndim = len(id_mm.variables['glamt'].shape)
        if ndim==2:
            kmaskt = id_mm.variables['tmask'][:,:]
            zlonF  = id_mm.variables['glamf'][:,:]
            zlatF  = id_mm.variables['gphif'][:,:]
            zlonT  = id_mm.variables['glamt'][:,:]
            zlatT  = id_mm.variables['gphit'][:,:]
            ze1T   = id_mm.variables['e1t'][:,:] / 1000. ; # km
            ze2T   = id_mm.variables['e2t'][:,:] / 1000. ; # km
            if l_CartesianCoordPresent:
                zYt = id_mm.variables[nm_y_t][:,:]
                zXt = id_mm.variables[nm_x_t][:,:]
                zYf = id_mm.variables[nm_y_f][:,:]
                zXf = id_mm.variables[nm_x_f][:,:]                
            if alsoF:
                kmaskf = id_mm.variables['fmask'][:,:]
            #
        elif ndim==4:
            kmaskt = id_mm.variables['tmask'][0,0,:,:]
            zlonF  = id_mm.variables['glamf'][0,:,:]
            zlatF  = id_mm.variables['gphif'][0,:,:]
            zlonT  = id_mm.variables['glamt'][0,:,:]
            zlatT  = id_mm.variables['gphit'][0,:,:]
            ze1T   = id_mm.variables['e1t'][0,:,:] / 1000. ; # km
            ze2T   = id_mm.variables['e2t'][0,:,:] / 1000. ; # km
            if l_CartesianCoordPresent:
                zYt = id_mm.variables[nm_y_t][0,:,:]
                zXt = id_mm.variables[nm_x_t][0,:,:]
                zYf = id_mm.variables[nm_y_f][0,:,:]
                zXf = id_mm.variables[nm_x_f][0,:,:]                
            if alsoF:
                kmaskf = id_mm.variables['fmask'][0,0,:,:]
        else:
            print(' ERROR [GetModelGrid]: unexpected number of dimmensions for variable `glamt` !')
            exit(0)
    #
    kmaskt  = np.array(kmaskt, dtype='i1')
    zlonT   = np.mod( zlonT, 360. )
    zlonF   = np.mod( zlonF, 360. )
    (nj,ni) = np.shape(kmaskt)
    #
    if not l_CartesianCoordPresent:
        # Conversion from Geographic coordinates (lat,lon) to Cartesian in km,
        #  ==> same North-Polar-Stereographic projection as RGPS data...        
        zXt = np.zeros((nj,ni))
        zYt = np.zeros((nj,ni))
        zXf = np.zeros((nj,ni))
        zYf = np.zeros((nj,ni))
        zYt[:,:], zXt[:,:] = ConvertGeo2CartesianNPSkm(zlatT, zlonT)
        zYf[:,:], zXf[:,:] = ConvertGeo2CartesianNPSkm(zlatF, zlonF)

    # Local resolution in km (for ):
    zResKM = np.zeros((nj,ni))
    zResKM[:,:] = np.sqrt( ze1T*ze1T + ze2T*ze2T )
    del ze1T, ze2T

    if alsoF:
        return kmaskt, zlatT, zlonT, zYt, zXt, zYf, zXf, zResKM, kmaskf, zlatF, zlonF
    else:
        return kmaskt, zlatT, zlonT, zYt, zXt, zYf, zXf, zResKM


def GetModelUVGrid( fNCmeshmask ):

    chck4f( fNCmeshmask)

    # Reading mesh metrics into mesh-mask file:
    with Dataset(fNCmeshmask) as id_mm:
        zlonV = id_mm.variables['glamv'][0,:,:]
        zlatV = id_mm.variables['gphiv'][0,:,:]
        zlonU = id_mm.variables['glamu'][0,:,:]
        zlatU = id_mm.variables['gphiu'][0,:,:]

    (nj,ni) = np.shape(zlonV)

    zXv = np.zeros((nj,ni))
    zYv = np.zeros((nj,ni))
    zXu = np.zeros((nj,ni))
    zYu = np.zeros((nj,ni))

    # Conversion from Geographic coordinates (lat,lon) to Cartesian in km,
    #  ==> same North-Polar-Stereographic projection as RGPS data...
    zlonV = np.mod( zlonV, 360. )
    zlonU = np.mod( zlonU, 360. )
    zYv[:,:], zXv[:,:] = ConvertGeo2CartesianNPSkm(zlatV, zlonV)
    zYu[:,:], zXu[:,:] = ConvertGeo2CartesianNPSkm(zlatU, zlonU)
    del zlatV, zlonV, zlatU, zlonU

    return zYv, zXv, zYu, zXu


def GetSeedMask( fFSmask, mvar='tmask' ):

    chck4f( fFSmask)

    # Reading mesh metrics into mesh-mask file:
    with Dataset(fFSmask) as id_mm:
        kmaskt = id_mm.variables[mvar][:,:]

    (nj,ni) = np.shape(kmaskt)

    kmaskt = np.array(kmaskt, dtype='i1')

    return kmaskt






def GetModelSeaIceConc( fNCsi3, name='siconc', krec=0, expected_shape=[] ):

    chck4f( fNCsi3)

    print('    * [GetModelSeaIceConc]: reading "'+name+'" at record '+str(krec)+' in '+fNCsi3+' !')
    with Dataset(fNCsi3) as id_si3:
        zsic   = id_si3.variables[name][krec,:,:]

    if len(expected_shape)>0:
        if np.shape(zsic) != expected_shape:
            print('ERROR [GetModelSeaIceConc]: wrong shape for sea-ice concentration read:',np.shape(zsic),', expected:',expected_shape)
            sys.exit(0)

    return zsic



def ncSaveCloudBuoys( cf_out, ptime, pIDs, pY, pX, pLat, pLon, mask=[], xtime=[],
                      tunits=tunits_default, fillVal=FillValue, corigin=None ):
    '''
    '''
    print('\n *** [ncSaveCloudBuoys]: About to generate file: '+cf_out+' ...')
    (Nt,) = np.shape(ptime)
    (Nb,) = np.shape(pIDs)
    if np.shape(pY)!=(Nt,Nb) or np.shape(pX)!=(Nt,Nb) or np.shape(pLat)!=(Nt,Nb) or np.shape(pLon)!=(Nt,Nb):
        print('ERROR [ncSaveCloudBuoys]: one of the 2D arrays has a wrong shape!!!')
        exit(0)
    lSaveMask = (np.shape(mask)  == (Nt,Nb))
    lSaveTime = (np.shape(xtime) == (Nt,Nb)) ; # time for each buoy!
    #
    f_out = Dataset(cf_out, 'w', format='NETCDF4')
    #
    # Dimensions:
    cd_time = 'time'
    cd_buoy = 'buoy'
    f_out.createDimension(cd_time, None)
    f_out.createDimension(cd_buoy, Nb  )
    #
    # Variables:
    v_time = f_out.createVariable(cd_time,     'i4',(cd_time,))
    v_buoy = f_out.createVariable(cd_buoy,     'i4',(cd_buoy,))
    v_bid  = f_out.createVariable('id_buoy',   'i8',(cd_buoy,))
    x_lat  = f_out.createVariable('latitude' , 'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    x_lon  = f_out.createVariable('longitude', 'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    x_ykm  = f_out.createVariable('y_pos' ,    'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    x_xkm  = f_out.createVariable('x_pos',     'f4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
    #
    v_time.units = tunits
    v_bid.units  = 'ID of buoy'
    x_lat.units = 'degrees north'
    x_lon.units = 'degrees south'
    x_ykm.units    = 'km'
    x_xkm.units    = 'km'
    #
    if lSaveMask:
        v_mask = f_out.createVariable('mask',  'i1',(cd_time,cd_buoy,),                      zlib=True, complevel=9)
    if lSaveTime:
        x_tim  = f_out.createVariable('time_pos', 'i4',(cd_time,cd_buoy,), fill_value=fillVal, zlib=True, complevel=9)
        x_tim.units = tunits
    #
    v_buoy[:] = np.arange(Nb,dtype='i8')
    v_bid[:]  = pIDs[:]
    #
    for jt in range(Nt):
        v_time[jt]  = ptime[jt]
        x_lat[jt,:] =  pLat[jt,:]
        x_lon[jt,:] =  pLon[jt,:]
        x_ykm[jt,:] =    pY[jt,:]
        x_xkm[jt,:] =    pX[jt,:]
        if lSaveMask:
            v_mask[jt,:] = mask[jt,:]
        if lSaveTime:
            x_tim[jt,:] = xtime[jt,:]
    #
    if corigin:
         f_out.Origin = corigin
    f_out.About  = 'Lagrangian sea-ice drift'
    f_out.Author = 'Generated with `'+path.basename(argv[0])+'` of `sitrack` (L. Brodeau, 2023)'
    f_out.close()
    print('      ===> '+cf_out+' saved!')
    #
    return 0

def LoadNCtime( cfile, ltime2d=False, iverbose=0 ):
    '''
       Open & inspect a `mojito-generated` NetCDF file (like that generated by `ncSaveCloudBuoys`)
       and load the time vector only !
    '''
    #
    chck4f(cfile)
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        if not 'time' in list_dim:
            print(' ERROR [LoadNCtime()]: no dimensions `'+time+'` found into input file!'); exit(0)
        #
        listv = list( id_in.variables.keys() )
        if not 'time' in listv:
            print(' ERROR [LoadNCtime()]: no variable `'+time+'` found into input file!'); exit(0)
        Nt = id_in.dimensions['time'].size

        # Time record:
        ctunits = id_in.variables['time'].units
        if not ctunits == tunits_default:
            print(' ERROR [LoadNCtime()]: we expect "'+tunits_default+'" as units for the time record vector, yet we have: '+ctunits); exit(0)
        print('    * [LoadNCtime] => reading "time" ('+str(Nt)+' records) in file '+path.basename(cfile))
        ztime = id_in.variables['time'][:]

        if ltime2d:
            if not 'time_pos' in listv:
                print(' ERROR [LoadNCtime()]: no variable `'+time_pos+'` found into input file!'); exit(0)
            ctunits = id_in.variables['time_pos'].units
            if not ctunits == tunits_default:
                print(' ERROR [LoadNCtime()]: we expect "'+tunits_default+'" as units for the 2D time record array, yet we have: '+ctunits); exit(0)
            print('    * [LoadNCtime] => reading "time_pos" in file '+path.basename(cfile))
            ztime2d = id_in.variables['time_pos'][:,:]
            (kt,_) = np.shape(ztime2d)
            if kt!=Nt:
                print(' ERROR [LoadNCtime()]: array `time_pos` has not the same number of records as `time`!!!'); exit(0)
            #
            return Nt, ztime, ztime2d
        else:
            return Nt, ztime

#time_pos:lolo

def LoadNCdata( cfile, krec=-1, lmask=False, lGetTimePos=False, iverbose=0 ):
    '''
       Open & inspect a `sitrack-generated` NetCDF file (like that generated by `ncSaveCloudBuoys`)
       and load the data

       * krec: record to extract, if krec==None => all records are extracted

    '''
    #
    listv_needed = ['id_buoy','latitude','longitude','y_pos','x_pos']
    if lGetTimePos:
        listv_needed.append('time_pos')

    chck4f(cfile)
    with Dataset(cfile) as id_in:

        list_dim = list( id_in.dimensions.keys() )
        for cd in ['time','buoy']:
            if not cd in list_dim:
                print(' ERROR [LoadNCdata()]: no dimensions `'+cd+'` found into input file!'); exit(0)

        listv = list( id_in.variables.keys() )
        for cv in listv_needed:
            if not cv in listv:
                print(' ERROR [LoadNCdata()]: no variable `'+cv+'` found into input file!'); exit(0)

        Nt = id_in.dimensions['time'].size
        nP = id_in.dimensions['buoy'].size
        if iverbose>0:
            print('\n *** Total number of records in the file = ', Nt)
            print('  *** Total number of buoys in the file = ', nP)

        # Time record:
        ctunits = id_in.variables['time'].units
        if not ctunits == tunits_default:
            print(' ERROR [LoadNCdata()]: we expect "'+tunits_default+'" as units for the time record vector, yet we have: '+ctunits)
            exit(0)

        if krec>=0:
            idxR =  krec  ; # => `idxR` is an integer !
        else:
            # All records will be read:
            idxR = np.arange( Nt, dtype=int ) ; # => `idxR` is a vector of integers !

        ztime = id_in.variables['time'][idxR]

        # Buoys' IDs:
        kBIDs    = np.zeros(nP, dtype=int)
        kBIDs[:] = id_in.variables['id_buoy'][:]

        # Coordinates:
        zlat  = id_in.variables['latitude'][idxR,:]
        zlon  = id_in.variables['longitude'][idxR,:]
        zy    = id_in.variables['y_pos'][idxR,:]
        zx    = id_in.variables['x_pos'][idxR,:]
        if lmask:
            zmsk = id_in.variables['mask'][idxR,:]
        if lGetTimePos:
            ztpos = id_in.variables['time_pos'][idxR,:]

    zlon[:] = np.mod(zlon, 360.) ; # Longitudes in the [0:360] frame...

    if krec>=0:
        zLatLon = np.zeros((nP,2))
        zYX     = np.zeros((nP,2))
        zLatLon[:,:] = np.array([zlat,zlon]).T
        zYX[:,:]     = np.array([ zy , zx ]).T
    else:
        zLatLon = np.zeros((Nt,nP,2))
        zYX     = np.zeros((Nt,nP,2))
        for jt in range(Nt):
            zLatLon[jt,:,:] = np.array([zlat[jt,:],zlon[jt,:]]).T
            zYX[jt,:,:]     = np.array([  zy[jt,:] , zx[jt,:]]).T

    if lmask:
        if lGetTimePos:
            return ztime, kBIDs, zLatLon, zYX, zmsk, ztpos
        else:
            return ztime, kBIDs, zLatLon, zYX, zmsk
    else:
        if lGetTimePos:
            return ztime, kBIDs, zLatLon, zYX, ztpos
        else:
            return ztime, kBIDs, zLatLon, zYX


def SeedFileTimeInfo( fSeedNc, ltime2d=False, iverbose=0 ):
    from re import split
    from math import ceil, floor
    #
    cSeed = str.replace( path.basename(fSeedNc), 'SELECTION_', '' )
    cSeed = str.replace( cSeed, '.nc', '' )
    cBtch = split('_',path.basename(fSeedNc))[2]
    #
    chck4f(fSeedNc)
    #
    if ltime2d:
        ntr, zt, zt2d  = LoadNCtime( fSeedNc, ltime2d=True, iverbose=iverbose )
        idate0, idateN = np.min(zt2d), np.max(zt2d)
    else:
        ntr, zt       = LoadNCtime( fSeedNc,               iverbose=iverbose )
        idate0, idateN = zt[0], zt[ntr-1]
        zt2d = []
    #
    if ntr>1:
        print('    * [SeedFileTimeInfo] => earliest and latest time position in the SEED file: '+e2c(idate0)+' - '+e2c(idateN))
        # Generous rounding at the hour sharp:
        idate0, idateN = int( floor(idate0/3600.)*3600. ), int( ceil(idateN/3600.)*3600. )
        print('    * [SeedFileTimeInfo]  ==> will actually use rounded to the hour! => '+e2c(idate0)+' - '+e2c(idateN))
    elif ntr==1:
        print('    * [SeedFileTimeInfo] => only 1 time record in the SEED file! time = ', e2c(idate0))
        idateN = idate0 ; # just to be sure
    else:
        print('ERROR [SeedFileTimeInfo]: problem!'); exit(0)
    #
    return idate0, idateN, cSeed, cBtch, zt2d



def ModelFileTimeInfo( fModelNc, iverbose=0 ):
    #
    from re import split
    #
    with Dataset(fModelNc) as ds_mod:
        cv_t = GetNameTimeDim( ds_mod ) ; # find the name for time/calendar vector:
        #
        Nt = ds_mod.dimensions[cv_t].size
        if ds_mod.variables[cv_t].units != tunits_default:
            print('    * [ModelFileTimeInfo] => "'+ds_mod.variables[cv_t].units+'" is the units for time in file',fModelNc)
            print('                          ==> need to convert to '+tunits_default+' !')
            time0 = ds_mod.variables[cv_t]
            ztime = ConvertTimeToEpoch(time0[:], time0.units, time0.calendar )
        else:
            ztime = np.array( ds_mod.variables[cv_t][:] , dtype='i4' )
    #
    print('    * [ModelFileTimeInfo] => '+str(Nt)+' records in input MODEL file!')
    #
    idate0, idateN = np.min(ztime), np.max(ztime)
    print('    * [ModelFileTimeInfo] => earliest and latest time position in the MODEL file: '+e2c(idate0)+' - '+e2c(idateN))
    #
    # Infer name of NEMO CONFIG and experiment from SI3 file:
    vn = split('_',path.basename(fModelNc))
    zz = split('-',vn[1])
    nconf = vn[0]
    if len(zz)==1:
        zz = split('-',vn[0])
        nconf = zz[0]
    nexpr = zz[1]
    print('    * [ModelFileTimeInfo] => Strings for `config` and `experiment` =>', nconf, nexpr, '\n')
    #    
    #print(' Nt, ztime, idate0, idateN, nconf, nexpr =', Nt, ztime, idate0, idateN, nconf, nexpr); exit(0)
    #
    return Nt, ztime, idate0, idateN, nconf, nexpr

