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
iplot=1

lUseActualTime = True ; # => if set to True, each buoy will be tracked for the exact same amount of time as its RGPS counterpart
#                       #    => otherwize, it is tracked for 3 days...

rdt = 3600 ; # time step (must be that of model output ice velocities used)

toDegrees = 180./pi

isubsamp_fig = 72 ; # frequency, in number of model records, we spawn a figure on the map (if idebug>2!!!)

iUVstrategy = 1 ; #  What U,V should we use inside a given T-cell of the model?
#                 #  * 0 => use the same MEAN velocity in the whole cell => U = 0.5*(U[j,i-1] + U[j,i]), V = 0.5*(V[j-1,i] + U[j,i])
#                 #  * 1 => use the same NEAREST velocity in the whole cell => U = U[@ nearest U-point], V = V[@ nearest V-point]

def __argument_parsing__():
    '''
    ARGUMENT PARSING / USAGE
    '''
    import argparse as ap
    #
    parser = ap.ArgumentParser(description='SITRACK ICE PARTICULES TRACKER')
    rqrdNam = parser.add_argument_group('required arguments')
    rqrdNam.add_argument('-i', '--fsi3', required=True,  help='output file of SI3 containing ice velocities ans co')
    rqrdNam.add_argument('-m', '--fmmm', required=True, help='model `mesh_mask` file of NEMO config used in SI3 run')
    rqrdNam.add_argument('-s', '--fsdg', required=True,  help='seeding file')
    #
    parser.add_argument('-k', '--krec' , type=int, default=0, help='record of seeding file to use to seed from')
    parser.add_argument('-e', '--dend' ,  default=None,       help='date at which to stop')
    #
    args = parser.parse_args()
    print('')
    print(' *** SI3 file to get ice velocities from => ', args.fsi3)
    print(' *** SI3 `mesh_mask` metrics file        => ', args.fmmm)
    print(' *** Seeding file and record to use      => ', args.fsdg, args.krec )
    if args.dend:
        print(' *** Overidding date at which to stop =>', args.dend )
    #
    return args.fsi3, args.fmmm, args.fsdg, args.krec, args.dend



if __name__ == '__main__':


    print('')
    print('##########################################################')
    print('#            SITRACK ICE PARTICULES TRACKER              #')
    print('##########################################################\n')


    cf_uv, cf_mm, fNCseed, jrecSeed, cdate_stop = __argument_parsing__()

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

    
    # Are we in a idealized seeding or not:
    fncSsplt = split('\.', fNCseedBN)[0]
    fncSsplt = split('_',fncSsplt)
    if fncSsplt[2] in ['nemoTsi3','nemoTmm']:
        print('\n *** Seems to be an idealized seeding of type "'+fncSsplt[2]+'"')
        cdtbin = '_idlSeed'
        csfkm = '_Xkm'
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

        
    creskm = csfkm[1:]
    print(' *** The spatial resolution and dt_bin infered from file name: '+creskm+', '+cdtbin)
    
    # Some strings and start/end date of Seeding input file:

    idateSeedA, idateSeedB, SeedName, SeedBatch, zTpos = sit.SeedFileTimeInfo( fNCseed, ltime2d=lUseActualTime, iverbose=idebug )
    
    if lUseActualTime:
        (n2,nB) = np.shape(zTpos)
        if n2!=2:
            print('ERROR: wrong shape for the 2D time array!'); exit(0)
        print(' ==> ok, we got the 2 time positions for '+str(nB)+' buoys! (`lUseActualTime==True`)\n')
        
    # Same for model input file + time records info:
    Nt0, ztime_model, idateModA, idateModB, ModConf, ModExp = sit.ModelFileTimeInfo( cf_uv, iverbose=idebug )
    
    # What records of model data can we use, based on time info from 2 input files above:
    date_stop = None
    if ldateStop:
        # A stop date explicitely asked at command line:
        #print('LOLO: cdate_stop =',cdate_stop)
        if len(cdate_stop)==19:
            date_stop = mjt.clock2epoch( cdate_stop )
        else:
            date_stop = mjt.clock2epoch( cdate_stop, precision='D', cfrmt='basic' )
        
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
    
    cfdir = './figs/tracking/'+creskm    
    if iplot>0 and not path.exists(cfdir):
        makedirs( cfdir, exist_ok=True )
    for cd in ['nc', 'npz' ]:
        if not path.exists(cd): mkdir(cd)

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
    cf_npz_itm = './npz/Initialized_buoys_'+SeedName+'.npz'


    ############################
    # Initialization / Seeding #
    ############################

    if path.exists(cf_npz_itm):
        # We save a lot of energy by using the previously generated intermediate backup file:        
        print('\n *** We found file '+cf_npz_itm+' here! So using it and skeeping first stage!')
        with np.load(cf_npz_itm) as data:
            nP    = data['nP']
            xPosG0 = data['xPosG0']
            xPosC0 = data['xPosC0']
            IDs   = data['IDs']
            vJIt  = data['vJIt']
            VRTCS = data['VRTCS']

    else:

        # Going through whole initialization / seeding process
        # ----------------------------------------------------

        with Dataset(cf_uv) as ds_UVmod:
            #xIC[:,:] = ds_UVmod.variables['siconc'][0,:,:] ; # We need ice conc. at t=0 so we can cancel buoys accordingly
            xIC[:,:] = ds_UVmod.variables['siconc'][kstrt,:,:] ; # We need ice conc. at `t=kstrt` so we can cancel buoys accordingly

        zt, zIDs, XseedG, XseedC = sit.LoadNCdata( fNCseed, krec=jrecSeed, iverbose=idebug )
        print('     => data used for seeding is read at date =',e2c(zt),'\n        (shape of XseedG =',np.shape(XseedG),')')
        
        (nP,_) = np.shape(XseedG)

        # We want an ID for each seeded buoy:
        IDs = np.array( range(nP), dtype=int) + 1 ; # Default! No ID=0 !!!


        IDs[:] = zIDs[:]
        del zIDs

        # Find the location of each seeded buoy onto the model grid:
        nP, xPosG0, xPosC0, IDs, vJIt, VRTCS = sit.SeedInit( IDs, XseedG, XseedC, xlatT, xlonT, xYf, xXf,
                                                             xResKM, imaskt, xIceConc=xIC, iverbose=idebug )
        del XseedG, XseedC

        # This first stage is fairly costly so saving the info:
        print('\n *** Saving intermediate data into '+cf_npz_itm+'!')
        np.savez_compressed( cf_npz_itm, nP=nP, xPosG0=xPosG0, xPosC0=xPosC0, IDs=IDs, vJIt=vJIt, VRTCS=VRTCS )

    del xResKM

    
    # ========= Should be an external function ===================================================
    # What is the first and last model record to use for each buoy:
    z1stModelRec = np.zeros(nP, dtype=int) + kstrt ; # default 1st record is used
    zLstModelRec = np.zeros(nP, dtype=int) + kstop ; # default 1st record is used

    if lUseActualTime:
        if nP != nB:
            print('ERROR: `nP != nB` !'); exit(0)                
        # Now, we need to check if any buoy initial time is actually beyond first record of the model:
        lLater = (zTpos[0,:]>=iTmA+int(rdt/2))
        if np.any(lLater):
            #print('WARNING: there are initial buoy time position that are beyond first model record:')
            (idxLate,) = np.where(lLater)
            for jb in idxLate:
                #print('    * buoy at pos.',jb,'starts at',e2c(zTpos[0,jb]))
                (idx,) = np.where(ztime_model+int(rdt/2)<zTpos[0,jb])
                jrc = idx[-1]+1
                #print(' for jrc =',jrc, '(instead of jrc =',kstrt,') we have ztime_model[jrc] =',e2c(ztime_model[jrc]))
                z1stModelRec[jb] = jrc

        # Same, but for last needed record:
        lEarlr = (zTpos[1,:]<iTmB-int(rdt/2))
        if np.any(lEarlr):
            print('WARNING: there are final buoy time position that are before last model record:')
            (idxEarl,) = np.where(lEarlr)
            for jb in idxEarl:
                print(' * buoy at pos.',jb,'ends at',e2c(zTpos[1,jb]),'(last used model time=',e2c(iTmB),')')
                (idx,) = np.where(ztime_model-int(rdt/2)>zTpos[1,jb])
                jrc = idx[0]-1
                print('  =>for jrc=',jrc, '(instead of jrc =',kstop,') we have ztime_model[jrc] =',e2c(ztime_model[jrc]))
                zLstModelRec[jb] = jrc
        #print('min and max for `zLstModelRec` =',np.min(zLstModelRec),np.max(zLstModelRec))
        #(idxMx,)=np.where(zLstModelRec==np.max(zLstModelRec))
        #print('Number of buoys that goes to the max:',len(idxMx))

        # DEBUG prompt! That's pretty important:
        for jb in range(0,nP,10):
            print(' * B '+str(jb)+': init. & end time=',e2c(zTpos[0,jb]),e2c(zTpos[1,jb]),', 1st & last model rec to use:',
                  e2c(ztime_model[z1stModelRec[jb]]),e2c(ztime_model[zLstModelRec[jb]]))

    # lili: ok we know which 1st record to use but not how to implement it in the loop to come yet! #lolo
    #exit(0)
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
    
    if lUseActualTime:
        k0 = z1stModelRec[:] - kstrt ; # It's a vector in this case!!!!  #fixme: it's only okay when dt_model=1h ???
    else:
        k0 = 0
    xPosC[k0,:,:] = xPosC0
    xPosG[k0,:,:] = xPosG0
    xmask[k0,:,:] = 1
    exit(0)
    del xPosC0, xPosG0

    if iplot>0 and idebug>1:
        mjt.ShowBuoysMap( 0, xPosG[0,:,1], xPosG[0,:,0], cfig=cfdir+'/INIT_Pos_buoys_'+SeedBatch+'_'+ModExp+'_'+'%4.4i'%(jt)+csfkm+'.png',
                          cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1., title='IceTracker: Init Seeding' ) ; #, pvIDs=IDs






    ######################################
    # Loop along model data time records #
    ######################################

    # Open Input SI3 data file
    ds_UVmod = Dataset(cf_uv)

    for jt in range(Nt):

        jrec = jt + kstrt ; # access into netCDF file...

        rtmod = ds_UVmod.variables['time_counter'][jrec] ; # time of model data (center of the average period which should = rdt)
        itime = int(rtmod - rdt/2.) ; # velocitie is average under the whole rdt, at the center!
        ctime = e2c(itime)
        print('\n *** Reading record #'+str(jrec+1)+'/'+str(Nt0)+' in SI3 file ==> date =',
              ctime,'(model:'+e2c(int(rtmod))+')')
        vTime[jt] = itime

        xIC[:,:] = ds_UVmod.variables['siconc'][jrec,:,:]
        xUu[:,:] = ds_UVmod.variables['u_ice'][jrec,:,:]
        xVv[:,:] = ds_UVmod.variables['v_ice'][jrec,:,:]

        print('   *   current number of buoys alive = '+str(iAlive.sum()))

        for jP in range(nP):

            if iAlive[jP]==1:

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

            ### if iAlive[jP]==1

        ### for jP in range(nP)

        # Updating in terms of lon,lat for all the buoys at once:
        xPosG[jt+1,:,:] = sit.CartNPSkm2Geo1D( xPosC[jt+1,:,:] )

        print('\n')
    ### for jt in range(Nt)

    ds_UVmod.close()

    vTime[Nt] = vTime[Nt-1] + int(rdt)


    # Masking arrays:
    xPosG = np.ma.masked_where( xmask==0, xPosG )
    xPosC = np.ma.masked_where( xmask==0, xPosC )

    # ==> time to save itime, xPosXX, xPosYY, xPosLo, xPosLa into a netCDF file !
    cdt1, cdt2 = split(':',e2c(vTime[0]))[0] , split(':',e2c(vTime[Nt]))[0] ; # keeps at the hour precision...
    cdt1, cdt2 = str.replace( cdt1, '-', '') , str.replace( cdt2, '-', '')
    cdt1, cdt2 = str.replace( cdt1, '_', 'h') , str.replace( cdt2, '_', 'h')
    corgn = 'NEMO-SI3_'+ModConf+'_'+ModExp
    cf_nc_out = './nc/'+corgn+'_tracking_'+SeedBatch+cdtbin+'_'+cdt1+'_'+cdt2+csfkm+'.nc'

    kk = sit.ncSaveCloudBuoys( cf_nc_out, vTime, IDs, xPosC[:,:,0], xPosC[:,:,1], xPosG[:,:,0], xPosG[:,:,1],
                               mask=xmask[:,:,0], corigin=corgn )



    if iplot>0:
        # Show on the map of the Arctic:
        for jt in range(Nt+1):
            if jt%isubsamp_fig == 0:
                zLon = np.ma.masked_where( xmask[jt,:,1]==0, xPosG[jt,:,1] )
                zLat = np.ma.masked_where( xmask[jt,:,0]==0, xPosG[jt,:,0] )
                ctag = cdt1+'-'+cdt2+'_'+'%4.4i'%(jt)
                cfig = cfdir+'/Pos_buoys_'+SeedBatch+csfkm+'_'+ModExp+'_'+ctag+'.png'

                mjt.ShowBuoysMap( vTime[jt], zLon, zLat, cfig=cfig,
                                  cnmfig=None, ms=5, ralpha=0.5, lShowDate=True, zoom=1.,
                                  title='IceTracker + SI3 '+ModExp+' u,v fields' ) ; # , pvIDs=IDs
                del zLon, zLat






    print('        => first and final dates in simulated trajectories:',e2c(vTime[0]),e2c(vTime[-1]),'\n')

