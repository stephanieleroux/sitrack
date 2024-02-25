#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:

'''

from sys import argv, exit
from os import path, makedirs, environ
import numpy as np
from re import split

import mojito   as mjt
import sitrack  as sit

#import random
from random import random, choices


idebug=0
iplot=0

#lRandomize = True ; rDssAmpl = 0.001

seeding_type='debug' ; # By default, ha

ldo_coastal_clean = False
distMax=100 ; distMin = 100 ; # how far from the nearest coast should our buoys be? [km]
fdist2coast_nc = 'dist2coast/dist2coast_4deg_North.nc'



l_CentralArctic = False ; # only keep points of the central Arctic

lAddFpoints = False


def __argument_parsing__():
    '''
    ARGUMENT PARSING / USAGE
    '''
    import argparse as ap
    #
    parser = ap.ArgumentParser(description='SITRACK ICE PARTICULES TRACKER')
    rqrdNam = parser.add_argument_group('required arguments')
    #
    rqrdNam.add_argument('-d', '--dat0', required=True,       help='initial date in the form <YYYY-MM-DD_hh:mm:ss>')
    #
    parser.add_argument('-m', '--fmmm' , default=None,        help='model `mesh_mask` file of NEMO config used in SI3 run')
    parser.add_argument('-i', '--fsi3' , default=None,        help='output file of SI3 containing sea-ice concentration')
    parser.add_argument('-v', '--nsic' , default='siconc',    help='name of sea-ice concentration in SI3 file (default="siconc")')
    parser.add_argument('-k', '--krec' , type=int, default=0, help='use sea-ice concentration at this record (only if you did specify a SI3 file with `--fsi3`!)')
    parser.add_argument('-S', '--ihss' , type=int, default=1, help='horizontal subsampling factor to apply')
    parser.add_argument('-f', '--fmsk' , default=None,        help='mask (on SI3 model domain) to control seeding region')
    parser.add_argument('-C', '--crsn' , type=int, default=0, help='apply this coarsening in km')
    parser.add_argument('-N', '--ncnf' , default='NANUK4',    help='name of the horizontak NEMO config used')
    parser.add_argument('--lsidfex' , type=int,default=0,     help='Switch to 1 for SIDFEX seeding.') 
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
    if args.crsn>=1:
        print(' *** Will apply a coarsening on cloud of points, at scale => ', args.crsn,'km' )
    if args.ncnf:
        print(' *** Name of the horizontak NEMO config used => ', args.ncnf)
    # SLX:
    if args.lsidfex:
        print(' *** SIDFEX seeding? => ', args.lsidfex)
    #
    return args.dat0, args.fsi3, args.nsic, args.krec, args.fmmm, args.ihss, args.fmsk, args.crsn, args.ncnf, args.lsidfex





if __name__ == '__main__':

    #xrand = []
    #for i in range(2000):
    #    rr = 2.*random()-1.
    #    xrand.append(rr)
    #    print(rr)

    #xrand = np.array(xrand)
    #print('* mean =', np.mean(xrand))
    #print('* max =', np.max(xrand))
    #print('* min =', np.min(xrand))
    #exit(0)


    
    print('')

    lnemoMM  = False
    lnemoSI3 = False
        
    cdate0, cf_si3, cv_sic, krec, cf_mm, iHSS, cf_force_msk, icrsn, CONF, lsidfex = __argument_parsing__()

    if cf_mm:
        lnemoMM  = True
        seeding_type = 'nemoTmm'
    if cf_si3:
        lnemoSI3 = True
        seeding_type = 'nemoTsi3'

    lForceSeedRegion = False
    if cf_force_msk:
        lForceSeedRegion = True

    lCoarsen = ( icrsn>=1 )
        
    if iplot>0:
        if   CONF=='NANUK4':
            name_proj = 'CentralArctic'
        elif CONF=='HUDSON4':
            name_proj = 'HudsonB'
        else:
            print('ERROR: CONF "'+CONF+'" is unknown (to know what proj to use for plots...)')
            exit(0)
    
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
    if lCoarsen:
        print('coarsening:',icrsn,'km')

    if iHSS<1 or iHSS>20:
        print('ERROR: chosen horizontal subsampling makes no sense iHSS=',iHSS)
        exit(0)        
    if lnemoSI3:
        if krec<0:
            print('ERROR: chosen record to read is < 0!',krec); exit(0)
    print('')

    if lCoarsen:
        # Coarsening:        
        if   icrsn==10:
            lAddFpoints=True
            #rd_ss =  7.5
            #rd_ss =  7. ; # to test !!!! => 10.12-10.14
            #rd_ss =  6.5 ; # to test !!!! => 10.08
            #rd_ss =  6.2 ; # to test !!!! => 10.05
            rd_ss =  6. ; # to test !!!! => fine...
            #zAmpRand = 0.05 ; # degrees
            zAmpRand = 0.02 ; # degrees
        elif icrsn==20:
            lAddFpoints=True
            rd_ss = 14.6 ; # real shit! 
            zAmpRand = 0.1 ; # degrees
        elif icrsn==40:
            rd_ss = 34.5 ; # real shit!
            zAmpRand = 0.15 ; # degrees
        elif icrsn==80:
            rd_ss = 74.75 ; # real shit!
            zAmpRand = 0.15 ; # degrees
        elif icrsn==160:
            rd_ss = 156. ; # real shit!
            zAmpRand = 0.15 ; # degrees
        elif icrsn==320:
            rd_ss = 315.6 ; # real shit!
            zAmpRand = 0.15 ; # degrees
        elif icrsn==640:
            rd_ss = 636. ; # real shit!
            zAmpRand = 0.15 ; # degrees
            ldo_coastal_clean=True; distMax=300 ; distMin = 100 ; # how far from the nearest coast should our buoys be? [km]
        else:
            print('ERROR: we do not know what `rd_ss` to pick for `icrsn` =',icrsn)
            exit(0)









    
            
    if lnemoMM or lnemoSI3:
        # Getting model grid metrics and friends:        
        if lAddFpoints:
            imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM, imaskF, xlatF, xlonF = sit.GetModelGrid( cf_mm, alsoF=True )
        else:
            imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM                       = sit.GetModelGrid( cf_mm )

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
        if lAddFpoints:
            XseedGC = sit.nemoSeed( imaskt, xlatT, xlonT, xIC, khss=iHSS, fmsk_rstrct=FSmask, platF=xlatF, plonF=xlonF )
        else:
            XseedGC = sit.nemoSeed( imaskt, xlatT, xlonT, xIC, khss=iHSS, fmsk_rstrct=FSmask )
        #
    elif seeding_type=='debug':
        if idebug in [0,1,2]: XseedGC = sit.debugSeeding()
        if idebug in [3]:     XseedGC = sit.debugSeeding1()
        #
     # SLX
    elif seeding_type=='sidfex':
        XseedGC,zIDs = sit.SidfexSeeding()
        print(XseedGC)
        print(zIDs)
        #
    else:
         print(' ERROR: `seeding_type` =',seeding_type,' is unknown!')
         exit(0)
         
    print('\n * Shape of XseedGC =',np.shape(XseedGC))
    
    (nP,_) = np.shape(XseedGC)



    # SLX   
    if not seeding_type=='sidfex':
          zIDs = np.array([ i+1 for i in range(nP) ] , dtype=int)

    zTime = np.array( [ mjt.clock2epoch(cdate0) ], dtype='i4' )    
    print('\n * Requested initialization date =', mjt.epoch2clock(zTime[0]))
    
    cdate = mjt.epoch2clock(zTime[0], precision='h')
    cdate = str.replace(cdate, '-', '')
    

    #if lRandomize:
    #    # Working with stereo projection rather than geographic coordinates to avoid North-Pole singularity:
    #    zYX = sit.Geo2CartNPSkm1D( XseedGC )
    #    #
    #    for jP in range(nP):
    #        zr1 = 2.*random()-1. ; # random number between -1 and 1
    #        zr2 = 2.*random()-1. ; # random number between -1 and 1
    #        #
    #        [zy,zx] = zYX[jP,:]
    #        #
    #        zy = zy + zAmpRand*110.*zr1 ; # because 1 degree ~ 110 km...
    #        zx = zx + zAmpRand*110.*zr2
    #        zYX[jP,:] = [zy,zx]
    #        #
    #    XseedGC[:,:] = sit.CartNPSkm2Geo1D( zYX )
    #    del zYX

        
    if ldo_coastal_clean:
        #
        cdata_dir = environ.get('DATA_DIR')
        if cdata_dir==None:
            print('\n ERROR: Set the `DATA_DIR` environement variable!\n'); exit(0)
        fdist2coast_nc = cdata_dir+'/data/dist2coast/dist2coast_4deg_North.nc'
        #

        #mjt.ShowBuoysMap( zTime[0], XseedGC[:,1], XseedGC[:,0],
        #                  cfig='01_BEFORE_LL_DIST2COAST.png', cnmfig=None, ms=5, ralpha=0.5, lShowDate=True,
        #                  zoom=1., title='Seeding initialization' )
        
        #mjt.ShowBuoysMap( zTime[0], XseedGC[:,1], XseedGC[:,0],
        #                  cfig='02_AFTER_LL_DIST2COAST.png', cnmfig=None, ms=5, ralpha=0.5, lShowDate=True,
        #                  zoom=1., title='Seeding initialization' )
        #exit(0)

        
        #if lRandomize and distMin<distMax:
        #    zrc = 2.*random()-1. ; # random number between -1 and 1
        #    MinDistFromLand =    0.5*(distMin + distMax) + zrc*0.5*(distMax - distMin)            
        #else:
        MinDistFromLand = 0.5*(distMin + distMax)
        print(' *** `MinDistFromLand` to be used =',MinDistFromLand,'km')
                    
        mask = mjt.MaskCoastal( XseedGC, rMinDistLand=MinDistFromLand, fNCdist2coast=fdist2coast_nc, convArray='C' )
        print(' * Need to remove '+str(nP-np.sum(mask))+' points because too close to land! ('+str(MinDistFromLand)+'km)')
        nP = np.sum(mask) ; # new size once buoys too close to land removed...
        (idxKeep,) = np.where(mask==1)
        zIDs = zIDs[idxKeep]
        XseedGC =  np.array( [ np.squeeze(XseedGC[idxKeep,0]), np.squeeze(XseedGC[idxKeep,1]) ] ).T
        del idxKeep
        
    
    #if l_CentralArctic:
    #    (idxKeep,) = np.where( (XseedGC[:,0]>70.) & (np.mod(XseedGC[:,1],360.)>120.) & (np.mod(XseedGC[:,1],360.)<280.) | (XseedGC[:,0]>82.) )
    #    nP = len(idxKeep)
    #    zIDs = zIDs[idxKeep]
    #    XseedGC =  np.array( [ np.squeeze(XseedGC[idxKeep,0]), np.squeeze(XseedGC[idxKeep,1]) ] ).T
    #    del idxKeep


    # Convert geoo coordinates to projected polar proj:
    XseedYX = sit.Geo2CartNPSkm1D( XseedGC ) ; # same for seeded initial positions, XseedGC->XseedYX


    creskm = ''
    if lCoarsen:
        creskm = str(icrsn)+'km'

    cextra = ''
    
    if iHSS>1:
        cextra = '_HSS'+str(iHSS)
        
    if lCoarsen:
        cextra='_'+creskm

        if idebug>0 and iplot>0:
            fdir = './figs/coarsen/'+creskm
            makedirs( fdir, exist_ok=True )
            ffig = fdir+'/sitrack_seeding_'+seeding_type+'_'+cdate+cextra+'_beforeCRSN'+'.png'
            cxtr = str.replace(cextra, '_', ' ')            
            mjt.ShowBuoysMap( zTime[0], XseedGC[:,1], XseedGC[:,0],
                              cfig=ffig, cnmfig=None, ms=5, ralpha=0.5, lShowDate=True,
                              zoom=1., title='Seeding initialization'+cxtr )
            
        #MIND: both XseedGC and XseedYX are in C-array-indexing...
        
        #if lRandomize:
        #    zrand = 2.*random()-1. ; # random number between -1 and 1
        #    rd_ss = rd_ss + zrand*rDssAmpl
            
        print('\n *** Applying spatial sub-sampling with radius: '+str(round(rd_ss,2))+'km')
        nPss, zYXss, idxKeep = mjt.SubSampCloud( rd_ss, XseedYX[:,:], convArray='C' )
        zGCss = XseedGC[idxKeep,:].copy()
        print('    ==> nP, nPss =',nP, nPss )
        del XseedGC, XseedYX
        #
        nP = nPss
        XseedGC = zGCss.copy()
        XseedYX = zYXss.copy()
        zIDs = zIDs[idxKeep]
        del zGCss, zYXss, idxKeep
        
    
        
    XseedGC = np.reshape( XseedGC, (1,nP,2) )
    XseedYX = np.reshape( XseedYX, (1,nP,2) )
    

    print('shape XseedYX =',np.shape(XseedYX))
    print('shape XseedGC =',np.shape(XseedGC))


    print('idx 0 =',XseedGC[0,::20,0])
    print('idx 1 =',XseedGC[0,::20,1])

    
    makedirs( './nc', exist_ok=True )

    #if lRandomize:
    #    import string
    #    crand = ''.join(choices(string.ascii_letters + string.digits, k=8))
    #    cextra += '_'+crand
    
    foutnc = './nc/sitrack_seeding_'+seeding_type+'_'+cdate+cextra+'.nc'    
    print('\n *** Saving seeding file for date =',mjt.epoch2clock(zTime[0]),'\n   => into:',foutnc)

    
    kk = sit.ncSaveCloudBuoys( foutnc, zTime, zIDs, XseedYX[:,:,0], XseedYX[:,:,1], XseedGC[:,:,0], XseedGC[:,:,1],
                               corigin='idealized_seeding' )

    if iplot>0:
        fdir = './figs/SEEDING'
        makedirs( fdir, exist_ok=True )
        ffig = fdir+'/sitrack_seeding_'+seeding_type+'_'+cdate+cextra+'.png'


        cextra = str.replace(cextra, '_', ' ')
        
        mjt.ShowBuoysMap( zTime[0], XseedGC[0,:,1], XseedGC[0,:,0],
                          cfig=ffig, nmproj=name_proj, cnmfig=None, ms=5, ralpha=0.5, lShowDate=True,
                          zoom=1., title='Seeding initialization'+cextra )
    

