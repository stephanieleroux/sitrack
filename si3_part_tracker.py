#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:

   *

'''

from sys import argv, exit
from os import path, mkdir, makedirs
import numpy as np
from re import split
from netCDF4 import Dataset
from math import atan2,pi

#import gonzag as gzg
import mojito  as mjt
from mojito import epoch2clock as e2c
import sitrack as sit


idebug=0
iplot=0

#lUse2DTime = False ; # => if set to True, each buoy will be tracked for the exact same amount of time as its RGPS counterpart
#                       #    => otherwize, it is tracked for 3 days...

rdt = 3600. ; #FIXME!!! time step [s] (must be that of model output ice velocities used)

toDegrees = 180./pi

ifreq_plot = 24 ; # frequency, in terms of number of model records, we spawn a figure on the map (if idebug>2!!!)

iUVstrategy = 1 ; #  What U,V should we use inside a given T-cell of the model?
#iUVstrategy = 0 ; #  What U,V should we use inside a given T-cell of the model?
#                 #  * 0 => use the same MEAN velocity in the whole cell => U = 0.5*(U[j,i-1] + U[j,i]), V = 0.5*(V[j-1,i] + V[j,i])
#                 #  * 1 => use the same NEAREST velocity in the whole cell => U = U[@ nearest U-point], V = V[@ nearest V-point]

def __argument_parsing__():
    '''
    ARGUMENT PARSING / USAGE
    '''    
    import argparse as ap
    global lUse2DTime
    #
    parser = ap.ArgumentParser(description='SITRACK ICE PARTICULES TRACKER')
    rqrdNam = parser.add_argument_group('required arguments')
    rqrdNam.add_argument('-i', '--fsi3', required=True,  help='output file of SI3 containing ice velocities ans co')
    rqrdNam.add_argument('-m', '--fmmm', required=True, help='model `mesh_mask` file of NEMO config used in SI3 run')
    rqrdNam.add_argument('-s', '--fsdg', required=True,  help='seeding file')
    #
    parser.add_argument('-k', '--krec' , type=int, default=0,  help='record of seeding file to use to seed from')
    parser.add_argument('-e', '--dend' , default=None,         help='date at which to stop')
    parser.add_argument('-F', '--fxdt',  action="store_true",  help='fixed tracking time (1D time array)')
    parser.add_argument('-N', '--ncnf' , default='NANUK4',     help='name of the horizontak NEMO config used')
    parser.add_argument('-p', '--plot' , type=int, default=0,  help='how often, in terms of model records, we plot the positions on a map')
    #
    args = parser.parse_args()
    print('')
    print(' *** SI3 file to get ice velocities from => ', args.fsi3)
    print(' *** SI3 `mesh_mask` metrics file        => ', args.fmmm)
    print(' *** Seeding file and record to use      => ', args.fsdg, args.krec )
    if args.dend:
        print(' *** Overidding date at which to stop =>', args.dend )        
    if args.ncnf:
        print(' *** Name of the horizontak NEMO config used => ', args.ncnf)
    #
    lUse2DTime = not args.fxdt
    
    return args.fsi3, args.fmmm, args.fsdg, args.krec, args.dend, args.ncnf, args.plot



if __name__ == '__main__':


    print('')
    print('##########################################################')
    print('#            SITRACK ICE PARTICULES TRACKER              #')
    print('##########################################################\n')


    cf_uv, cf_mm, fNCseed, jrecSeed, cdate_stop, CONF, ifreq_plot = __argument_parsing__()

    print(' lUse2DTime =',lUse2DTime)
    print('cf_uv =',cf_uv)
    print('cf_mm =',cf_mm)
    print('fNCseed =',fNCseed)
    print('jrecSeed =',jrecSeed)
    ldateStop=False
    if cdate_stop:
        ldateStop=True
        print('cdate_stop =',cdate_stop)
    print('\n')

    fNCseedBN = path.basename(fNCseed)

    if iplot>0:
        if   CONF=='NANUK4':
            name_proj = 'CentralArctic'
        elif CONF=='HUDSON4':
            name_proj = 'HudsonB'
        else:
            print('ERROR: CONF "'+CONF+'" is unknown (to know what proj to use for plots...)')
            exit(0)

    lplot = (ifreq_plot>0)
    if lplot:
        print(' *** We shall plot a map with the position of virtual buoys every '+str(ifreq_plot)+' model records => '+str(int(rdt)*ifreq_plot/3600)+' hours')


    csfkm, creskm = '', ''
    # Are we in a idealized seeding or not:
    fncSsplt = split('\.', fNCseedBN)[0]
    fncSsplt = split('_',fncSsplt)
    print(fncSsplt[2])
    if fncSsplt[2] in ['nemoTsi3','nemoTmm','sidfex']:
        print('\n *** Seems to be an idealized seeding of type "'+fncSsplt[2]+'"')
        cdtbin = '_idlSeed'        
        for ii in [1,2,3]:
            ckm = '_'+fncSsplt[-ii]
            if ckm[-2:]=='km':
                csfkm = ckm
                break
        #
    else:
        # Info about resolution from the seeding file name?
        lok=False
        itst=1        
        while not lok:
            itst-=1
            csfkm = '_'+split( '_', split('\.',fNCseedBN)[-2] )[itst]
            kres=-3+itst            
            cdtbin = '_'+split( '_', split('\.',fNCseedBN)[-2] )[kres]
            lok = ( csfkm[-2:]=='km' and cdtbin[1:3]=='dt' ) or ( cdtbin[1:3]=='dt' and itst==0 )
            if itst<-4:
                print('ERROR: we could not figure out `csfkm` and `cdtbin` from file name!',csfkm, cdtbin); exit(0)
        if itst==0: csfkm = ''


    if csfkm != '':
        creskm = csfkm[1:]
        print(' *** The spatial resolution and dt_bin infered from file name: '+creskm+', '+cdtbin)
    else:
        creskm = ''
    
    # Some strings and start/end date of Seeding input file:
    idateSeedA, idateSeedB, SeedName, SeedBatch, zTpos = sit.SeedFileTimeInfo( fNCseed, ltime2d=lUse2DTime, iverbose=idebug )
    print("====== SLX")
    print(zTpos) 


   # Same for model input file + time records info:
    Nt0, ztime_model, idateModA, idateModB, ModConf, ModExp = sit.ModelFileTimeInfo( cf_uv, iverbose=idebug )
    
    # What records of model data can we use, based on time info from 2 input files above:
    date_stop = None
    if ldateStop:
        # A stop date explicitely asked at command line:
        if len(cdate_stop)==19:
            date_stop = mjt.clock2epoch( cdate_stop )
        else:
            date_stop = mjt.clock2epoch( cdate_stop, precision='D', cfrmt='guess' )
        
    elif idateSeedB - idateSeedA >= 3600.:
        # => there are more than 1 record in the file we use for seeding!!!
        #    ==> which means we are likely to do replicate the exact same thing using the model data
        #    ==> so we stop at `idateSeedB` !!!
        date_stop = idateSeedB
    
        
    #
    Nt, kstrt, kstop, iTmA, iTmB = sit.GetTimeSpan( rdt, ztime_model, idateSeedA, idateModA, idateModB , iStop=date_stop )
    #
    if Nt<1:
        print(' QUITTING since no matching model records!')
        exit(0)

    crk = creskm
    cc = split('-',creskm)
    if len(cc)==2: crk = cc[1]
    cfdir = './figs/tracking/'+crk    
    if iplot>0 and not path.exists(cfdir):
        makedirs( cfdir, exist_ok=True )
    for cd in [ 'seed', 'nc', 'npz' ]:
        makedirs( cd, exist_ok=True )

    # Getting model grid metrics and friends:
    imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM = sit.GetModelGrid( cf_mm )

    if iUVstrategy==1:
        # Get extra U,V-point metrics:
        xYv, xXv, xYu, xXu = sit.GetModelUVGrid( cf_mm )
        
    # Allocating arrays for model data:
    (Nj,Ni) = np.shape( imaskt )
    xUu = np.zeros((Nj,Ni))
    xVv = np.zeros((Nj,Ni))
    xIC = np.zeros((Nj,Ni)) ; # Sea-ice concentration

    # We need a name for the intermediate backup file:
    cf_npz_itm = './seed/Initialized_buoys_'+SeedName+'_'+CONF+'.npz'

    ############################
    # Initialization / Seeding #
    ############################

    if path.exists(cf_npz_itm):
        # We save a lot of energy by using the previously generated intermediate backup file:        
        print('\n *** We found file '+cf_npz_itm+' here! So using it and skipping first stage!')
        with np.load(cf_npz_itm) as data:
            nP    = data['nP']
            xPosG0 = data['xPosG0']
            xPosC0 = data['xPosC0']
            IDs   = data['IDs']
            vJIt  = data['vJIt']
            VRTCS = data['VRTCS']
            idxK  = data['idxKeep']
            
    else:
        print('\n *** We did not find file '+cf_npz_itm+' ! => going through nearest point scanning...')
        # Going through whole initialization / seeding process
        # ----------------------------------------------------

        with Dataset(cf_uv) as ds_UVmod:
            xIC[:,:] = ds_UVmod.variables['siconc'][kstrt,:,:] ; # We need ice conc. at `t=kstrt` so we can cancel buoys accordingly

        zt, zIDs, XseedG, XseedC = sit.LoadNCdata( fNCseed, krec=jrecSeed, iverbose=idebug )
        print('     => data used for seeding is read at date =',e2c(zt),'\n        (shape of XseedG =',np.shape(XseedG),')')
        
        (nP,_) = np.shape(XseedG)

        # We want an ID for each seeded buoy:
        IDs = np.array( range(nP), dtype=int) + 1 ; # Default! No ID=0 !!!

        IDs[:] = zIDs[:]
        del zIDs

        # Find the location of each seeded buoy onto the model grid:
        nPn, xPosG0, xPosC0, IDs, vJIt, VRTCS, idxK = sit.SeedInit( IDs, XseedG, XseedC, xlatT, xlonT, xYf, xXf,
                                                                    xResKM, imaskt, xIceConc=xIC, iverbose=idebug )
        del XseedG, XseedC

        if nPn<nP:
            print('\n *** `SeedInit()` had to cancel '+str(nP-nPn)+
                  ' buoys! => adjusting `zTpos` and updating nP from '+str(nP)+' to '+str(nPn)+'!')
            #zTpos = zTpos[:,idxK]
            nP = nPn
        
        # This first stage is fairly costly, so saving the info:
        print('\n *** Saving intermediate data into '+cf_npz_itm+'!')
        np.savez_compressed( cf_npz_itm, nP=nP, xPosG0=xPosG0, xPosC0=xPosC0, IDs=IDs, vJIt=vJIt, VRTCS=VRTCS, idxKeep=idxK )

    ### if path.exists(cf_npz_itm)
    
    del xResKM

    
    # ========= Should be an external function ===================================================
    # What is the first and last model record to use for each buoy:
    z1stModelRec = np.zeros(nP, dtype=int) + kstrt ; # default 1st record is used
    zLstModelRec = np.zeros(nP, dtype=int) + kstop ; # default 1st record is used

    if lUse2DTime:

        (n2,nB) = np.shape(zTpos)
        if n2!=2 or nP!=nB:
            print('ERROR: wrong shape for the 2D time array `zTpos`! `n2,nB`, vs `nP`:',n2,nB,nP)
            exit(0)
        print(' ==> ok, we got the 2 time positions for '+str(nB)+' buoys! (`lUse2DTime==True`)\n')
        del nB
        
        zstarting_dates = np.unique( zTpos[0,:] )
        zending_dates   = np.unique( zTpos[1,:] )
        print('All different starting dates found:')
        for kd in zstarting_dates:
            (idx,) = np.where(zTpos[0,:]==kd)
            print('    *',e2c(kd),' ==> ',len(idx),'occurences!')
        print('\nAll different ending dates found:')
        for kd in zending_dates:
            (idx,) = np.where(zTpos[1,:]==kd)
            print('    *',e2c(kd),' ==> ',len(idx),'occurences!')
        print('')
        
        # Now, we need to check if any buoy initial time is actually beyond first record of the model:
        lLater = (zTpos[0,:]>=iTmA+int(rdt/2))
        if np.any(lLater):
            print('WARNING: there are initial buoy time position that are after first model record:')
            (idxLate,) = np.where(lLater)
            for jb in idxLate:
                print(' * buoy at pos.',jb,'stars at',e2c(zTpos[0,jb]),'(last global model time=',e2c(iTmB),')')
                #print('    * buoy at pos.',jb,'starts at',e2c(zTpos[0,jb]))
                (idx,) = np.where(ztime_model+int(rdt/2)<zTpos[0,jb])
                jrc = idx[-1]+1
                print('  => for jrc =',jrc, '(instead of jrc =',kstrt,') we have ztime_model[jrc] =',e2c(ztime_model[jrc]))
                z1stModelRec[jb] = jrc

        print('')
        # Same, but for last needed record:
        lEarlr = (zTpos[1,:]<iTmB-int(rdt/2))
        if np.any(lEarlr):
            print('WARNING: there are final buoy time position that are before last model record:')
            (idxEarl,) = np.where(lEarlr)
            for jb in idxEarl:
                print(' * buoy at pos.',jb,'ends at',e2c(zTpos[1,jb]),'(last global model time=',e2c(iTmB),')')
                (idx,) = np.where(ztime_model-int(rdt/2)>zTpos[1,jb])
                jrc = idx[0]-1
                print('  =>for jrc=',jrc, '(instead of jrc =',kstop,') we have ztime_model[jrc] =',e2c(ztime_model[jrc]))
                zLstModelRec[jb] = jrc

        if idebug>1:
            # Debug summary:
            for jb in range(0,nP,10):
                print(' * B '+str(jb)+': init. & end time=',e2c(zTpos[0,jb]),e2c(zTpos[1,jb]),', 1st & last model rec to use:',
                      e2c(ztime_model[z1stModelRec[jb]]),e2c(ztime_model[zLstModelRec[jb]]))

    #==============================================================================================================

    
    # Allocation for nP buoys:
    iAlive = np.zeros(      nP , dtype='i1') + 1 ; # tells if a buoy is alive (1) or zombie (0) (discontinued)
    vTime  = np.zeros( Nt+1, dtype=int ) ; # UNIX epoch time associated to position below
    xmask  = np.zeros((Nt+1,nP,2), dtype='i1')
    xPosC  = np.zeros((Nt+1,nP,2)) + sit.FillValue  ; # x-position of buoy along the Nt records [km]
    xPosG  = np.zeros((Nt+1,nP,2)) + sit.FillValue
    vMesh  = np.zeros((nP,4,2))   ; # stores for each buoy the 4 points defining the cell (coordinates of 4 surrounding F-points)
    lStillIn = np.zeros(nP, dtype=bool) ; # tells if a buoy is still within expected mesh/cell..

    # Initial values for some arrays:
    if lUse2DTime:
        xTime  = np.zeros((Nt+1,nP),dtype=int) + sit.FillValue ; # will contain the time of the model record used!!!
        for jb in range(nP):
            k0 = z1stModelRec[jb] - kstrt ; #fixme: it's only okay when dt_model=1h ???
            xPosC[k0,jb,:] = xPosC0[jb,:]
            xPosG[k0,jb,:] = xPosG0[jb,:]
            xTime[k0,jb]   = ztime_model[z1stModelRec[jb]] - int(rdt/2)
            xmask[k0,jb,:] = 1
    else:
        xPosC[0,:,:] = xPosC0
        xPosG[0,:,:] = xPosG0
        xmask[0,:,:] = 1

    if iplot>0 and idebug>0:
        mjt.ShowBuoysMap( 0, xPosG0[:,1], xPosG0[:,0], cfig=cfdir+'/INIT_Pos_buoys_'+SeedBatch+'_'+ModExp+csfkm+'.png',
                          nmproj=name_proj, cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1., title='IceTracker: Init Seeding' ) ; #, pvIDs=IDs

    del xPosC0, xPosG0



    ######################################
    # Loop along model data time records #
    ######################################

    # Open Input SI3 data file
    ds_UVmod = Dataset(cf_uv)

    for jt in range(Nt):

        jrec = jt + kstrt ; # access into netCDF file...

        itmod = int( ds_UVmod.variables['time_counter'][jrec] ); # time of model data (center of the average period which should = rdt)
        itime = itmod - int(rdt/2.) ; # velocitie is average under the whole rdt, at the center!
        ctime = e2c(itime)
        print('\n *** Reading record #'+str(jrec+1)+'/'+str(Nt0)+' in SI3 file ==> date =',
              ctime,'(model:'+e2c(itmod)+')')
        vTime[jt] = itime

        xIC[:,:] = ds_UVmod.variables['siconc'][jrec,:,:]
        xUu[:,:] = ds_UVmod.variables['u_ice'][jrec,:,:]
        xVv[:,:] = ds_UVmod.variables['v_ice'][jrec,:,:]

        print('   *   current number of buoys alive = '+str(iAlive.sum()))

        for jP in range(nP):            
            
            if iAlive[jP]==1 and jrec>=z1stModelRec[jP] and jrec<=zLstModelRec[jP]:

                [ ry  , rx   ] = xPosC[jt,jP,:] ; # km !
                [ rlat, rlon ] = xPosG[jt,jP,:] ; # degrees !

                [jnT,inT] = vJIt[jP,:]; # indices for T-point at center of current cell

                if idebug>0:
                    print('\n    * BUOY ID:'+str(IDs[jP])+' => jt='+str(jt)+': ry, rx =', ry, rx,'km'+': rlat, rlon =', rlat, rlon )

                ######################### N E W   M E S H   R E L O C A T I O N #################################
                if not lStillIn[jP]:
                    if idebug>0: print('      +++ RELOCATION NEEDED for buoy with ID:'+str(IDs[jP])+' +++')

                    [ [ jbl, jbr, jur, jul ], [ ibl, ibr, iur, iul ] ] = VRTCS[jP,:,:]

                    if idebug>0:
                        print('     ==> 4 corner points of our mesh (anti-clockwise, starting from BLC) =',
                              [ [jbl,ibl],[jbr,ibr],[jur,iur],[jul,iul] ])

                    # The current mesh/cell:
                    vMesh[jP,:,:] = [ [xYf[jbl,ibl],xXf[jbl,ibl]], [xYf[jbr,ibr],xXf[jbr,ibr]],
                                      [xYf[jur,iur],xXf[jur,iur]], [xYf[jul,iul],xXf[jul,iul]] ]

                ### if not lStillIn[jP]
                ###########################################################################################################################

                #if idebug>2:
                #    # We can have a look in the mesh:
                #    cnames    = np.array([ 'P'+str(i+1)+': '+str(VRTCS[jP,0,i])+','+str(VRTCS[jP,1,i]) for i in range(4) ], dtype='U32')
                #    #sit.PlotMesh( (rlat,rlon), xlatF, xlonF, VRTCS[jP,:,:].T, vnames=cnames,
                #    #              fig_name='mesh_lon-lat_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                #    #              pcoor_extra=(xlatT[jnT,inT],xlonT[jnT,inT]), label_extra='T-point' )
                #    gzg.PlotMesh( ( ry , rx ),  xYf ,  xXf,  VRTCS[jP,:,:].T, vnames=cnames,
                #                  fig_name='mesh_X-Y_buoy'+'%3.3i'%(jP)+'_jt'+'%4.4i'%(jt)+'.png',
                #                  pcoor_extra=(xYt[jnT,inT],xXt[jnT,inT]), label_extra='T-point' )


                # j,i indices of the cell we are dealing with = that of the F-point aka the upper-right point !!!
                [jT,iT] = vJIt[jP,:]

                # ASSUMING THAT THE ENTIRE CELL IS MOVING AT THE SAME VELOCITY: THAT OF U-POINT OF CELL
                # zU, zV = xUu[jT,iT], xVv[jT,iT] ; # because the F-point is the upper-right corner
                if   iUVstrategy == 0:
                    zU = 0.5*(xUu[jT,iT]+xUu[jT,iT-1])
                    zV = 0.5*(xVv[jT,iT]+xVv[jT-1,iT])
                    #
                elif iUVstrategy == 1:
                    # If the segment that goes from our buoy position to the F-point of the cell
                    # intersects the segment that joins the 2 V-points of the cell (v[j-1,i],v[j,i]),
                    # then it means that the nearest U-point is the one at `i-1` !
                    Fpnt = [xYf[jT,iT],xXf[jT,iT]] ; # y,x coordinates of the F-point of the cell
                    llum1 = sit.intersect2Seg( [ry,rx], Fpnt,  [xYv[jT-1,iT],xXv[jT-1,iT]], [xYv[jT,iT],xXv[jT,iT]] )
                    llvm1 = sit.intersect2Seg( [ry,rx], Fpnt,  [xYu[jT,iT-1],xXu[jT,iT-1]], [xYu[jT,iT],xXu[jT,iT]] )
                    if llum1:
                        zU = xUu[jT,iT-1]
                    else:
                        zU = xUu[jT,iT]
                    if llvm1:
                        zV = xVv[jT-1,iT]
                    else:
                        zV = xVv[jT,iT]
                    if idebug>1:
                        print( ' ++ Buoy position is:',ry,rx)
                        print( ' ++ position of lhs & rhs U-point:',xYu[jT,iT-1],xXu[jT,iT-1], xYu[jT,iT],xXu[jT,iT], ' llum1=',llum1)
                        print( ' ++ position of lower & upper V-point:',xYv[jT,iT-1],xXv[jT,iT-1], xYv[jT,iT],xXv[jT,iT], ' llvm1=',llvm1)

                if idebug>0:
                    print('    =>> read velocity at ji,jj=',iT,jT)
                    print('    * ice velocity of the mesh: u,v =',zU, zV, 'm/s')

                # Displacement during the upcomming time step:
                dx = zU*rdt
                dy = zV*rdt
                if idebug>0: print('      ==> displacement during `dt`: dx,dy =',dx,dy, 'm')

                # => position [km] of buoy after this time step will be:
                rx_nxt = rx + dx/1000. ; # [km]
                ry_nxt = ry + dy/1000. ; # [km]
                xPosC[jt+1,jP,:] = [ ry_nxt, rx_nxt ]
                xmask[jt+1,jP,:] = [    1  ,    1   ]
                
                if lUse2DTime:
                    xTime[jt+1,jP] = itime + rdt ; # at `jt+1`, time is itime+rdt ! (xTime[0,:] filled earlier)
                
                # Is it still inside our mesh:
                lSI = sit.IsInsideQuadrangle( ry_nxt, rx_nxt, vMesh[jP,:,:] )
                lStillIn[jP] = lSI
                if idebug>0: print('      ==> Still inside the same mesh???',lSI)

                if not lSI:
                    # => point is exiting the current cell!

                    # Tells which of the 4 cell walls the point has crossed:
                    icross = sit.CrossedEdge( [ry,rx], [ry_nxt,rx_nxt], VRTCS[jP,:,:], xYf, xXf, iverbose=idebug )

                    # Tells in which adjacent cell the point has moved:
                    inhc = sit.NewHostCell( icross, [ry,rx], [ry_nxt,rx_nxt], VRTCS[jP,:,:], xYf, xXf,  iverbose=idebug )

                    # Update the mesh indices according to the new host cell:
                    VRTCS[jP,:,:],vJIt[jP,:] = sit.UpdtInd4NewCell( inhc, VRTCS[jP,:,:], vJIt[jP,:] )

                    # Based on updated new indices, some buoys might get killed:
                    icncl = sit.Survive( IDs[jP], vJIt[jP,:], imaskt, pIceC=xIC,  iverbose=idebug )
                    if icncl>0: iAlive[jP]=0

                ### if not lSI

            ### if iAlive[jP]==1 and jrec>=z1stModelRec[jP] and jrec<=zLstModelRec[jP]

        ### for jP in range(nP)

        # Updating in terms of lon,lat for all the buoys at once:
        xPosG[jt+1,:,:] = sit.CartNPSkm2Geo1D( xPosC[jt+1,:,:] )

        print('\n')
    ### for jt in range(Nt)

    ds_UVmod.close()

    vTime[Nt] = vTime[Nt-1] + int(rdt)


    # Masking arrays:
    #xPosG = np.ma.masked_where( xmask==0, xPosG )
    #xPosC = np.ma.masked_where( xmask==0, xPosC )
    #xTime = np.ma.masked_where( xmask[:,:,0]==0, xTime )

    
    # ==> time to save itime, xPosXX, xPosYY, xPosLo, xPosLa into a netCDF file !
    cdt1, cdt2 = split(':',e2c(vTime[0]))[0] , split(':',e2c(vTime[Nt]))[0] ; # keeps at the hour precision...
    cdt1, cdt2 = str.replace( cdt1, '-', '') , str.replace( cdt2, '-', '')
    cdt1, cdt2 = str.replace( cdt1, '_', 'h') , str.replace( cdt2, '_', 'h')
    corgn = 'NEMO-SI3_'+ModConf+'_'+ModExp

    if not lUse2DTime:
        # Save series at each model time step:
        cf_nc_out = './nc/'+corgn+'_tracking_'+SeedBatch+cdtbin+'_'+cdt1+'_'+cdt2+csfkm+'.nc'
        kk = sit.ncSaveCloudBuoys( cf_nc_out, vTime, IDs, xPosC[:,:,0], xPosC[:,:,1], xPosG[:,:,0], xPosG[:,:,1],
                                   mask=xmask[:,:,0], corigin=corgn )



    # Now we should create the 2-record (initial and final) nc file (mandatory if `lUse2DTime` !):

    z2XY, z2GC, zMSK = np.zeros((2,nP,2)), np.zeros((2,nP,2)), np.zeros((2,nP,2),dtype='i1') ; # [record,n.buoys,yx]    
    
    if lUse2DTime:
        
        if idebug>1:
            for jb in range(0,nP,10):
                print(' * xPosC at jt=Nt-3,Nt-2,Nt-1,Nt =',
                      Nt-3,Nt-2,Nt-1,Nt,'  ',xPosC[Nt-3,jb,0],xPosC[Nt-2,jb,0],xPosC[Nt-1,jb,0],xPosC[Nt,jb,0],'zLstModelRec[jb]-kstrt=',zLstModelRec[jb]-kstrt)
            print()
        
        zTim = np.zeros((2,nP),dtype=int)
        for jb in range(nP):
            # First valid record for buoy jb:
            k0 = z1stModelRec[jb] - kstrt ; #fixme: it's only okay when dt_model=1h ???
            z2XY[0,jb,:] = xPosC[k0,jb,:]
            z2GC[0,jb,:] = xPosG[k0,jb,:]
            zTim[0,jb]   = xTime[k0,jb]
            zMSK[0,jb,:] = xmask[k0,jb,:]
            # Last valid record for buoy jb:
            kN = zLstModelRec[jb] - kstrt + 1; # yes! `+1` is needed! #fixme: it's only okay when dt_model=1h ???
            z2XY[1,jb,:] = xPosC[kN,jb,:]
            z2GC[1,jb,:] = xPosG[kN,jb,:]
            zTim[1,jb]   = xTime[kN,jb]
            zMSK[1,jb,:] = xmask[kN,jb,:]
            #
        zvt = np.array([ np.mean(zTim[0,:]), np.mean(zTim[1,:]) ])
        
    else:
        z2XY[0,:,:] = xPosC[0,:,:]
        z2GC[0,:,:] = xPosG[0,:,:]
        zMSK[0,:,:] = xmask[0,:,:]
        z2XY[1,:,:] = xPosC[Nt,:,:]
        z2GC[1,:,:] = xPosG[Nt,:,:]
        zMSK[1,:,:] = xmask[Nt,:,:]
        zTim = [] ;  # not needed!
        #
        zvt = np.array([ vTime[0], vTime[Nt] ])

    
    # Save the netCDF files with 2 records (first and last) with exact same pattern as done with RGPS data in mojito:
    cdt1, cdt2 = split(':',e2c(zvt[0]))[0] , split(':',e2c(zvt[1]))[0] ; # keeps at the hour precision...
    cdt1, cdt2 = str.replace( cdt1, '-', '') , str.replace( cdt2, '-', '')
    cdt1, cdt2 = str.replace( cdt1, '_', 'h') , str.replace( cdt2, '_', 'h')    
    cf_nc_out = './nc/'+corgn+'_tracking12_'+SeedBatch+cdtbin+'_'+cdt1+'_'+cdt2+csfkm+'.nc'
    
    kk = sit.ncSaveCloudBuoys( cf_nc_out, zvt, IDs, z2XY[:,:,0], z2XY[:,:,1], z2GC[:,:,0], z2GC[:,:,1],
                               mask=zMSK[:,:,0], xtime=zTim, corigin=corgn )
    
    if iplot>0:
        # Show first and last valid records on the map of the Arctic:
        for jt in range(2):
            zLon = np.ma.masked_where( zMSK[jt,:,1]==0, z2GC[jt,:,1] )
            zLat = np.ma.masked_where( zMSK[jt,:,0]==0, z2GC[jt,:,0] )
            ctag = cdt1+'-'+cdt2+'_'+'%4.4i'%(jt)
            cfig = cfdir+'/Pos_buoys_1stLst_'+SeedBatch+csfkm+'_'+ModExp+'_'+ctag+'.png'

            mjt.ShowBuoysMap( zvt[jt], zLon, zLat, cfig=cfig, nmproj=name_proj, 
                              cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1.,
                              title='IceTracker + SI3 '+ModExp+' u,v fields' ) ; # , pvIDs=IDs
            del zLon, zLat
                

    if lplot:
        print('\n *** Will now generate the maps!')
        # Show on the map of the Arctic:
        for jt in range(Nt+1):
            if jt%ifreq_plot == 0:
                zLon = np.ma.masked_where( xmask[jt,:,1]==0, xPosG[jt,:,1] )
                zLat = np.ma.masked_where( xmask[jt,:,0]==0, xPosG[jt,:,0] )
                ctag = cdt1+'-'+cdt2+'_'+'%4.4i'%(jt)
                cfig = cfdir+'/Pos_buoys_'+SeedBatch+csfkm+'_'+ModExp+'_'+ctag+'.png'

                mjt.ShowBuoysMap( vTime[jt], zLon, zLat, cfig=cfig, nmproj=name_proj,
                                  cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1.,
                                  title='IceTracker + SI3 '+ModExp+' u,v fields' ) ; # , pvIDs=IDs
                del zLon, zLat




    print('        => global first and final dates in simulated trajectories:',e2c(zvt[0]),e2c(zvt[1]),'\n')

