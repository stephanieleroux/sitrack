#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:

'''

from sys import argv, exit
from os import path, mkdir, makedirs, environ
import numpy as np
from re import split

import mojito   as mjt
import sitrack  as sit

idebug=0
iplot=1

seeding_type='debug' ; # By default, ha

iHSS=15 ; # horizontal sub-sampling for NEMO points initialization...

ldo_coastal_clean = False
MinDistFromLand  = 200 ; # how far from the nearest coast should our buoys be? [km]
fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'

l_CentralArctic = False ; # only keep points of the central Arctic







def __argument_parsing__():
    '''
    ARGUMENT PARSING / USAGE
    '''
    import argparse as ap
    #
    parser = ap.ArgumentParser(description='SITRACK ICE PARTICULES TRACKER')
    rqrdNam = parser.add_argument_group('required arguments')
    #
    rqrdNam.add_argument('-d', '--dat0', required=True,  help='initial date in the form <YYYY-MM-DD_hh:mm:ss>')
    #
    parser.add_argument('-m', '--fmmm' , default=None,   help='model `mesh_mask` file of NEMO config used in SI3 run')
    parser.add_argument('-i', '--fsi3' , default=None,   help='output file of SI3 containing sea-ice concentration')
    parser.add_argument('-v', '--nsic' , default='siconc',help='name of sea-ice concentration in SI3 file (default="siconc")')
    parser.add_argument('-k', '--krec' , type=int, default=0,   help='record of seeding file to use to seed from (only if use SI3 file!)')
    parser.add_argument('-S', '--ihss' , type=int, default=1,   help='horizontal subsampling factor to apply')
    parser.add_argument('-f', '--fmsk' , default=None,   help='mask (on SI3 model domain) to control seeding region')
    args = parser.parse_args()

    if args.fsi3 and not args.fmmm:
        #print('ERROR: chose between SI3 output file or MeshMask file !!! (i.e. `-m` or `-i`)')
        print('ERROR: you have to specify a MeshMask file with `-m` when using SI3 file!')
        exit(0)
    
    print('')
    print(' *** Date for initilization => ', args.dat0)
    if args.fsi3:
        print(' *** SI3 file to get sea-ice concentration from => ', args.fsi3)
        print('     ==> name of sea-ice concentration field    => ', args.nsic )
        print('     ==> record to use                   => ', args.krec )
    if args.fmmm:
        print(' *** SI3 `mesh_mask` metrics file        => ', args.fmmm)
    print(' *** Horizontal subsampling factor to apply => ', args.ihss)
    if args.fmsk:
        print(' *** Will apply masking on seeding data, file to use => ', args.fmsk )
    #
    return args.dat0, args.fsi3, args.nsic, args.krec, args.fmmm, args.ihss, args.fmsk





if __name__ == '__main__':

    print('')
    if ldo_coastal_clean:
        cdata_dir = environ.get('DATA_DIR')
        if cdata_dir==None:
            print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
        fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'

    lnemoMM  = False
    lnemoSI3 = False
        
    cdate0, cf_si3, cv_sic, krec, cf_mm, iHSS, cf_force_msk = __argument_parsing__()

    if cf_mm:
        lnemoMM  = True
        seeding_type = 'nemoTmm'
    if cf_si3:
        lnemoSI3 = True
        seeding_type = 'nemoTsi3'

    lForceSeedRegion = False
    if cf_force_msk:
        lForceSeedRegion = True

    
    print('\n *** Input SUMMARY:')
    print('cdate0 =',cdate0)
    if lnemoSI3:
        print('cf_si3 =',cf_si3)
        print('cv_sic =',cv_sic)
        print('krec =',krec)
    if lnemoMM:
        print('cf_mm =',cf_mm)
    print('iHSS =',iHSS)
    if lForceSeedRegion:
        print('cf_force_msk =',cf_force_msk)
    if iHSS<1 or iHSS>20:
        print('ERROR: chosen horizontal subsampling makes no sense iHSS=',iHSS)
        exit(0)        
    if lnemoSI3:
        if krec<0:
            print('ERROR: chosen record to read is < 0!',krec); exit(0)
    print('')

            
    if lnemoMM or lnemoSI3:
        # Getting model grid metrics and friends:
        imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM = sit.GetModelGrid( cf_mm )

    if lnemoSI3:
        xIC = sit.GetModelSeaIceConc( cf_si3, name=cv_sic, krec=krec, expected_shape=np.shape(imaskt) )
        #
    elif lnemoMM:
        # A fake sea-ice concentration:
        xIC = np.ones( np.shape(imaskt) )

    FSmask = None
    if lForceSeedRegion:
        FSmask = sit.GetSeedMask( cf_force_msk, mvar='tmask' )
        print(FSmask[::50,::20])                                                                                                                                                               
        if np.shape(FSmask) != np.shape(imaskt):
            print('ERROR: `shape(FSmask) != shape(imaskt)`'); exit(0)

    
    ############################
    # Initialization / Seeding #
    ############################


    if seeding_type in ['nemoTmm','nemoTsi3']:
        XseedG = sit.nemoSeed( imaskt, xlatT, xlonT, xIC, khss=iHSS, fmsk_rstrct=FSmask )
        #
    elif seeding_type=='debug':
        if idebug in [0,1,2]: XseedG = sit.debugSeeding()
        if idebug in [3]:     XseedG = sit.debugSeeding1()
        #
    else:
         print(' ERROR: `seeding_type` =',seeding_type,' is unknown!')
         exit(0)
         
    print('\n * Shape of XseedG =',np.shape(XseedG))
    
    (nP,_) = np.shape(XseedG)




    zIDs = np.array([ i+1 for i in range(nP) ] , dtype=int)

    zTime = np.array( [ mjt.clock2epoch(cdate0) ], dtype='i4' )    
    print('\n * Requested initialization date =', mjt.epoch2clock(zTime[0]))
    
    print( zIDs  )
    print( zTime )

    if ldo_coastal_clean:
        mask = mjt.MaskCoastal( XseedG, rMinDistLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc, convArray='C' )
        print(' * Need to remove '+str(nP-np.sum(mask))+' points because too close to land! ('+str(MinDistFromLand)+'km)')
        nP = np.sum(mask) ; # new size once buoys too close to land removed...
        (idxKeep,) = np.where(mask==1)
        zIDs = zIDs[idxKeep]
        XseedG =  np.array( [ np.squeeze(XseedG[idxKeep,0]), np.squeeze(XseedG[idxKeep,1]) ] ).T
        del idxKeep
        
    
    if l_CentralArctic:
        #zYXkm = mjt.Geo2CartNPSkm1D( XseedG, lat0=50., lon0=180. )        
        #zdistP = np.sqrt( zYXkm[:,0]*zYXkm[:,0] + zYXkm[:,1]*zYXkm[:,1] )
        #print(zdistP[::10])
        #(idxKeep,) = np.where(zdistP<1000.)
        #(idxKeep,) = np.where( (XseedG[:,0]>82.) | ((np.mod(XseedG[:,1],360.)>120.) & (np.mod(XseedG[:,1],360.)<280.)) )
        (idxKeep,) = np.where( (XseedG[:,0]>70.) & (np.mod(XseedG[:,1],360.)>120.) & (np.mod(XseedG[:,1],360.)<280.) | (XseedG[:,0]>82.) )
        nP = len(idxKeep)
        zIDs = zIDs[idxKeep]
        XseedG =  np.array( [ np.squeeze(XseedG[idxKeep,0]), np.squeeze(XseedG[idxKeep,1]) ] ).T
        del idxKeep

        
    XseedC = sit.Geo2CartNPSkm1D( XseedG ) ; # same for seeded initial positions, XseedG->XseedC
        
    XseedG = np.reshape( XseedG, (1,nP,2) )
    XseedC = np.reshape( XseedC, (1,nP,2) )
    
    makedirs( './nc', exist_ok=True )
    foutnc = './nc/sitrack_seeding_'+seeding_type+'_'+mjt.epoch2clock(zTime[0], precision='D')+'_HSS'+str(iHSS)+'.nc'
    
    print('\n *** Saving seeding file for date =',mjt.epoch2clock(zTime[0]))
    
    kk = sit.ncSaveCloudBuoys( foutnc, zTime, zIDs, XseedC[:,:,0], XseedC[:,:,1], XseedG[:,:,0], XseedG[:,:,1],
                               corigin='idealized_seeding' )

    if iplot>0:

        makedirs( './figs', exist_ok=True )
        ffig = './figs/sitrack_seeding_'+seeding_type+'_'+mjt.epoch2clock(zTime[0], precision='D')+'.png'

        
        mjt.ShowBuoysMap( 0, XseedG[0,:,1], XseedG[0,:,0],
                          cfig=ffig, cnmfig=None, ms=5, ralpha=0.5, lShowDate=True,
                          zoom=1., title='Seeding initialization' )
        # pvIDs=zIDs,
    

