
import numpy as np

rmin_conc = 0.1 ; # ice concentration below which we disregard the point...
rFoundKM = 2.5


def GetTimeSpan( dt, vtime_mod, iSdA, iMdA, iMdB, iStop=None, iverbose=0 ):
    '''
         * iStop: date at which to stop !
                  => if not provided we stop at latest model date!!!
    '''
    ltStop = False
    if iStop: ltStop = True

    # Confrontation of time info seeding and model:
    from .util import epoch2clock as e2c
    #
    if iSdA < iMdA-dt/2 or iSdA> iMdB-dt/2:
        print('PROBLEM: time in the seeding file ('+e2c(iSdA)+') is outside of what model spans!')
        exit(0)
    #
    kt0 = np.argmin(np.abs(vtime_mod[:]-iSdA))
    if iSdA >= vtime_mod[kt0]: kt0 += 1
    itM0 = vtime_mod[kt0]
    print('    * [GetTimeSpan]: First record needed =',kt0,'of model input file =>',e2c(itM0))
    #
    if ltStop:
        ktN = np.argmin(np.abs(vtime_mod[:]-iStop))
    else:
        ktN = len(vtime_mod) - 1
    #
    itMN = vtime_mod[ktN]
    print('    * [GetTimeSpan]: Last record needed =',ktN,'of model input file =>',e2c(itMN))
    Nt = ktN - kt0 + 1
    print('       ==> '+str(Nt)+' model records')
    print('       ==> that makes '+str(round((itMN-itM0)/(3600*24),3))+' days of ice particule tracking.')
    #
    return Nt, kt0, ktN, itM0, itMN






def _ccw_( pcA, pcB, pcC ):
    '''
        * pcA: coordinates of point A => [y_A,x_A]
        * etc...
    '''
    return (pcC[0]-pcA[0])*(pcB[1]-pcA[1]) > (pcB[0]-pcA[0])*(pcC[1]-pcA[1])

def intersect2Seg( pcA, pcB, pcC, pcD ):
    '''
     Return true if line segments AB and CD intersect

        * pcA: coordinates of point A => [y_A,x_A]
        * etc...
    '''
    return ( _ccw_(pcA,pcC,pcD) != _ccw_(pcB,pcC,pcD) ) and ( _ccw_(pcA,pcB,pcC) != _ccw_(pcA,pcB,pcD) )



def Survive( kID, kjiT, pmskT, pIceC=[],  iverbose=0 ):
    '''
        Suite of tests to decide whether a buoy should be killed or not...
    '''
    (Nj,Ni) = np.shape(pmskT) ; # shape of the NEMO domain...
    
    [jT,iT] = kjiT
    ikill   = 0

    # Test: too close to NEMO domain boundaries:
    if ikill==0:
        if jT in [0,1,Nj-2,Nj-1] or iT in [0,1,Ni-2,Ni-1]:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (reaching NEMO domain boundaries)')
    
    # Test on land-sea mask / continent:
    if ikill==0:
        zmt = pmskT[jT,iT] + pmskT[jT,iT+1]+pmskT[jT+1,iT]+pmskT[jT,iT-1]+pmskT[jT-1,iT-1]
        if zmt < 5:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (too close to or over the land-sea mask)')

    # Test on sea-ice concentration:
    if ikill==0:
        if len(np.shape(pIceC))==2:
            zic = 0.2*(pIceC[jT,iT] + pIceC[jT,iT+1]+pIceC[jT+1,iT]+pIceC[jT,iT-1]+pIceC[jT-1,iT-1])
            if iverbose>1: print('     =>> 5P sea-ice concentration =',zic)
        if zic < rmin_conc:
            ikill = ikill+1
            if iverbose>0: print('        ===> I CANCEL buoy '+str(kID)+'!!! (mean 5P ice concentration at T-point of the cell =',zic,')')
    #
    return ikill




def SeedInit( pIDs, pSG, pSC, platT, plonT, pYf, pXf, pResolKM, maskT, xIceConc=[], iverbose=0 ):
    '''
    '''
    from time    import time
    from .locate import NearestPoint, FindContainingCell
    #
    (nP,n2) = np.shape(pSG)
    #
    if np.shape(pSC)!=(nP,n2):
        print('ERROR [SeedInit]: shape disagreement for `pSG` and `pSC`!'); exit(0)
    if n2!=2:
        print('ERROR [SeedInit]: wrong shape for `pSG` and `pSC`!'); exit(0)
    #
    #                  jF,iF:       # j,i indices of the F-point that defines the current mesh/cell
    #                                 # (was the nearest point when we searched for nearest point,
    #                                 #  as buoys move moves within the cell, it might not be the nearest point)
    zjiT    = np.zeros((nP,2),   dtype=int)
    zJIvrt  = np.zeros((nP,2,4), dtype=int)
    kmask = np.zeros( nP ,     dtype='i1') + 1

    time_start = time()
    
    for jP in range(nP):

        if jP%100==0:
            print('   * [SeedInit()]: finding nearest T-point and host cell... buoy #'+str(jP)+' / '+str(nP)+' ...')
        
        if iverbose>1:
            print('   * [SeedInit()]: focus on buoy with ID:'+str(pIDs[jP]))

        # Initial position of the buoy:
        zlat,zlon  = pSG[jP,0], pSG[jP,1] ; # degrees!
        zy  ,zx    = pSC[jP,0], pSC[jP,1] ; # km !
        zic = 1.
        
        # 1/ Nearest T-point on NEMO grid:
        (jT, iT) = NearestPoint( (zlat,zlon), platT, plonT, rd_found_km=rFoundKM, resolkm=pResolKM, max_itr=10 )

        if jT<0 or iT<0:
            kmask[jP] = 0
            if iverbose>0: print('        ===> I CANCEL buoy '+str(pIDs[jP])+'!!! (NO nearest T-point found for ',zlat,zlon,')')
                            
        if kmask[jP] == 1:            
            # Ok a nearest point was found!    
            if iverbose>1:
                print('     ==> nearest T-point for ',zlat,zlon,' on NEMO grid:', jT, iT, '==> lat,lon:',
                      round(platT[jT,iT],3), round(plonT[jT,iT],3))
            # Tests for canceling buoy or not:
            icncl = Survive( pIDs[jP], [jT,iT] , maskT, pIceC=xIceConc, iverbose=iverbose )
            if icncl>0:
                kmask[jP] = 0
                if iverbose>0: print('        ===> I CANCEL buoy '+str(pIDs[jP])+'!!! (Did not pass `Survive()` test!)')
            
        if kmask[jP] == 1:
            # Everything is okay, now we locate the cell/mesh (polygon joining 4 F-points) that includes
            # our target point (zy,zx)
            lPin, zjiT[jP,:], zJIvrt[jP,:,:] = FindContainingCell( (zy,zx), (jT,iT), pYf, pXf, iverbose=iverbose )
            #
            if not lPin:
                print('WARNING [SeedInit()]: could not find the proper F-point cell!!!')
                print('         => when lookin for point:',zlat,zlon)
                kmask[jP] = 0
                if iverbose>0: print('        ===> I CANCEL buoy '+str(pIDs[jP])+'!!! (NO proper F-point cell found)')
                        
    ### for jP in range(nP)
    time_stop = time()
    print(' * [SeedInit]: number of seconds it took to locate all the points in a target grid cell:', time_stop-time_start )
    
    # Now we can shrink the arrays based on `kmask`
    iKeep = np.arange(nP,dtype=int)
    nPn   = np.sum(kmask)
    if nPn<nP:
        (iGone,) = np.where(kmask==0)
        print(' * [SeedInit()]: '+str(nP-nPn)+' "to-be-seeded" buoys have to be canceled.')
        print('          => their IDs:',pIDs[iGone])
        print('        ==> need to shrink some arrays, number of valid buoys is now',nPn)
        nP = nPn
        (iKeep,) = np.where(kmask==1)
        
                
    return nP, pSG[iKeep,:], pSC[iKeep,:], pIDs[iKeep], zjiT[iKeep,:], zJIvrt[iKeep,:,:], iKeep



def CrossedEdge( pP1, pP2, ji4vert, pY, pX,  iverbose=0 ):
    '''
        * pP1    : point 1 as [y1,x1]
        * pP2    : point 2 as [y2,x2]
        * ji4vert : array(2,4), containg the 4 j,i indices of the 4 vertices of the mesh
    '''
    #
    for kk in range(4):
        kp1 = (kk+1)%4
        [j1,i1] = ji4vert[:,kk]
        [j2,i2] = ji4vert[:,kp1]
        #
        ll = intersect2Seg( pP1, pP2,  [pY[j1,i1],pX[j1,i1]], [pY[j2,i2],pX[j2,i2]] )
        if ll: break
        #
    if iverbose>0:
        vdir = ['bottom', 'right-hand', 'upper', 'left-hand']
        print('    [CrossedEdge()]: particle is crossing the '+vdir[kk]+' edge of the mesh!')
    return kk+1


def NewHostCell( kcross, pP1, pP2, ji4vert, pY, pX,  iverbose=0 ):
    '''
       Find in which adjacent cell, the point has moved into !

       * kcross : integer describing which (old) cell edge the point has crossed

       In rare cases, the buoy could possibly move into diagonally adjacent cells
       (not only bottom,right,upper,left cells...)
       These "diagonally-adjacent" meshes yields output: 5,6,7,8
    '''
    #
    [ [ jbl, jbr, jur, jul ], [ ibl, ibr, iur, iul ] ] = ji4vert[:,:]
    #
    knhc = kcross
    if  kcross==1:
        # Crosses the bottom edge:
        if   intersect2Seg( pP1, pP2, [pY[jbl,ibl],pX[jbl,ibl]], [pY[jbl-1,ibl],pX[jbl-1,ibl]] ):
            knhc=5 ; # bottom left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jbr,ibr],pX[jbr,ibr]], [pY[jbr-1,ibr],pX[jbr-1,ibr]] ):
            knhc=6 ; # bottom right diagonal
        #
    elif kcross==2:
        # Crosses the RHS edge:
        if   intersect2Seg( pP1, pP2, [pY[jbr,ibr],pX[jbr,ibr]], [pY[jbr,ibr+1],pX[jbr,ibr+1]] ):
            knhc=6 ; # bottom right diagonal
        elif intersect2Seg( pP1, pP2, [pY[jur,iur],pX[jur,iur]], [pY[jur,iur+1],pX[jur,iur+1]] ):
            knhc=7 ; # bottom right diagonal
        #
    elif kcross==3:
        # Crosses the upper edge:
        if   intersect2Seg( pP1, pP2, [pY[jul,iul],pX[jul,iul]], [pY[jul+1,iul],pX[jul+1,iul]] ):
            knhc=8 ; # upper left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jur,iur],pX[jur,iur]], [pY[jur+1,iur],pX[jur+1,iur]] ):
            knhc=7 ; # bottom right diagonal
        #
    elif kcross==4:
        # Crosses the LHS edge:
        if   intersect2Seg( pP1, pP2, [pY[jul,iul],pX[jul,iul]], [pY[jul,iul-1],pX[jul,iul-1]] ):
            knhc=8 ; # upper left diagonal
        elif intersect2Seg( pP1, pP2, [pY[jbl,ibl],pX[jbl,ibl]], [pY[jbl,ibl-1],pX[jbl,ibl-1]] ):
            knhc=5 ; # bottom left diagonal
    #
    if iverbose>0:
        vdir = ['bottom', 'RHS', 'upper', 'LHS', 'bottom-LHS', 'bottom-RHS', 'upper-RHS', 'upper-LHS' ]
        print('    *** Particle is moving into the '+vdir[knhc-1]+' mesh !')
    #
    return knhc



def UpdtInd4NewCell( knhc, ji4vert, kjiT, iverbose=0 ):
    '''
        Update the mesh indices according to the new host cell
    '''
    if   knhc==1:
        # went down:
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[0] = kjiT[0]-1
    elif knhc==5:
        if iverbose>0: print(' * [UpdtInd4NewCell()]: WE HAVE A 5 !!!!')
        # went left + down
        ji4vert[1,:] = ji4vert[1,:] - 1
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[1] = kjiT[1]-1
        kjiT[0] = kjiT[0]-1
    elif knhc==6:
        if iverbose>0: print(' * [UpdtInd4NewCell()]: WE HAVE A 6 !!!!')
        # went right + down
        ji4vert[1,:] = ji4vert[1,:] + 1
        ji4vert[0,:] = ji4vert[0,:] - 1
        kjiT[1] = kjiT[1]+1
        kjiT[0] = kjiT[0]-1
    elif knhc==2:
        # went to the right
        ji4vert[1,:] = ji4vert[1,:] + 1
        kjiT[1] = kjiT[1]+1
    elif knhc==3:
        # went up:
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[0] = kjiT[0]+1
    elif knhc==7:
        if iverbose>0: print(' * [UpdtInd4NewCell()]: WE HAVE A 7 !!!!')
        # went right + up
        ji4vert[1,:] = ji4vert[1,:] + 1
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[1] = kjiT[1]+1
        kjiT[0] = kjiT[0]+1
    elif knhc==8:
        if iverbose>0: print(' * [UpdtInd4NewCell()]: WE HAVE A 8 !!!!')
        # went right + up
        ji4vert[1,:] = ji4vert[1,:] - 1
        ji4vert[0,:] = ji4vert[0,:] + 1
        kjiT[1] = kjiT[1]-1
        kjiT[0] = kjiT[0]+1
    elif knhc==4:
        # went to the left
        ji4vert[1,:] = ji4vert[1,:] - 1
        kjiT[1] = kjiT[1]-1
    else:
        print('ERROR: unknown direction, knhc=',knhc)
        exit(0)

    return ji4vert, kjiT





# Misc. seeding functions (when seeding not based on input file)...

def debugSeeding():
    zLatLon = np.array([
        [84. ,  20.],
        [89. ,  50.],
        [63. , -11.], # point outside of domain
        [85. , 100.],        
        [89. , 100.],
        [79. , 180.],
        [76. ,  46.], # No ice, Barents Sea
        [75. , 190.],
        [85.2, -15.], # Too close to land
        [75. , 210.],
        [75. , -72.], # Baffin bay
        [83. , 200.],
        [79. , -42.], # over mask (Greenland!)
        [85. , 300.]  ])
    return zLatLon

def debugSeeding1():
    zLatLon = np.array([ [75.,190.] ])
    return zLatLon

def SidfexSeeding( filepath='./sidfexloc.dat' ):
  ''' Generate seeding lon lat from SIDFEX buoys read from text file.
  Returns an array with lat lon and another with buoy IDs.
  '''
  print("============= SIDFEX SEEDING")
  #debug example zLatLon = np.array([ [75.,190.] ])
  sidfexdat = ReadFromSidfexDatFile( filepath )
  print('sidfexdat', sidfexdat)
  print('sidfexdat.shape', sidfexdat.shape)
  print('sidfexdat.ndim', sidfexdat.ndim)
  if sidfexdat.ndim > 1:
      zLatLon = sidfexdat[:,[2,1]] # reverse order of columns so that it is now lat lon.
      zIDs = sidfexdat[:,0].astype(int)
  else:
      zLatLonsingle = sidfexdat[[2,1]]
      zLatLon=np.stack((zLatLonsingle,zLatLonsingle))  # here we duplicate the line to avoid problem when read by script generate_seeding
      zIDssingle = sidfexdat[0].astype(int)
      zIDs = np.stack((zIDssingle,zIDssingle)) # here we duplicate the line to avoid problem when read by script generate_seeding
  return zLatLon,zIDs
    


  def ReadFromSidfexDatFile( filepath='./sidfexloc.dat' ):
  '''
      Reads from text file sidfex.dat and returns an array with id, lon, lat for all SIDFEX buoys at a given date.
  '''
  from os.path import exists
  if (exists(filepath)):
      # reads buoys id and locations (lon lat) from text file
      listbuoys = np.genfromtxt(open(filepath))
      # Keep track of how many buoys are read in input file
      NBUOYTOTR=listbuoys.shape[0]
      print("======== Nb of buoys read from sidfex file......"+str(NBUOYTOTR))
  else:
      print("============================================")
      print("=== Error ! 'sidfexloc.dat' file is missing.")
      print("============================================")
      exit(0)
  return listbuoys

def nemoSeed( pmskT, platT, plonT, pIC, khss=1, fmsk_rstrct=[],
              platF=[], plonF=[] ):

    lAddF = ( np.shape(platF)==np.shape(pmskT) and np.shape(plonF)==np.shape(pmskT) )

    #if lAddF:
    #    # we create our custom F mask (1 where all 4 suround T points have 1!)
    #    pmskF = np.zeros(np.shape(pmskT), dtype='i1')
    #    pmskF[1:-1,1:-1] = ( pmskT[2:,1:-1] + pmskT[1:-1,2:] + pmskT[:-2,1:-1] + pmskT[1:-1,:-2] ) / 4
    
    zmsk = pmskT[::khss,::khss]
    zlat = platT[::khss,::khss]
    zlon = plonT[::khss,::khss]

    if lAddF:
        #zmskF = pmskF[::khss,::khss]
        zlatF = platF[::khss,::khss]
        zlonF = plonF[::khss,::khss]
        #del pmskF
        
    (Nj,Ni) = np.shape(zmsk)
    ztmp = np.zeros((Nj,Ni))

    lfmsk_rstrct = (np.shape(fmsk_rstrct) == np.shape(pmskT))
    if lfmsk_rstrct:
        maskR = fmsk_rstrct[::khss,::khss]
        if np.shape(maskR) != (Nj,Ni):
            print('ERROR [nemoSeed()]: restricted area mask does not agree in shape with model output!'); exit(0)
    
    msk_T = np.zeros((Nj,Ni), dtype='i1')
    msk_T[:,:] = zmsk[:,:]
    #if lAddF:
    #    msk_F = np.zeros((Nj,Ni), dtype='i1')
    #    msk_F[:,:] = zmskF[:,:]

    if lfmsk_rstrct:
        msk_T[:,:] = msk_T[:,:]*maskR[:,:]
        #if lAddF:
        #    msk_F[:,:] = msk_F[:,:]*maskR[:,:]
            
    # Only north of ...
    msk_T[np.where(zlat < 55.)] = 0
    #if lAddF:
    #    msk_F[np.where(zlat < 55.)] = 0
        
    # Only over a decent concentration of ice:
    ztmp[:,:] = pIC[::khss,::khss]
    msk_T[np.where(ztmp < 0.9)] = 0
    #if lAddF:
    #    msk_F[np.where(ztmp < 0.9)] = 0
            
    (idy_keep, idx_keep) = np.where( msk_T==1 )
    NbIPt = len(idy_keep)
    zLatLon = np.zeros((NbIPt,2))

    for jp in range(NbIPt):
        jj = idy_keep[jp]
        ji = idx_keep[jp]
        zLatLon[jp,:] = [ zlat[jj,ji],zlon[jj,ji] ]
        
    
    if lAddF:
        # we create our custom F mask (1 where all 4 suround T points have 1!)
        msk_F = np.zeros(np.shape(msk_T), dtype='i1')
        msk_F[1:-1,1:-1] = ( msk_T[2:,1:-1] + msk_T[1:-1,2:] + msk_T[:-2,1:-1] + msk_T[1:-1,:-2] ) / 4
        #
        (idy_keepF, idx_keepF) = np.where( msk_F==1 )
        NbIPtF = len(idy_keepF)
        zLatLonF = np.zeros((NbIPtF,2))

        for jp in range(NbIPtF):
            jj = idy_keepF[jp]
            ji = idx_keepF[jp]
            zLatLonF[jp,:] = [ zlatF[jj,ji],zlonF[jj,ji] ]
        print(' * [nemoSeed()]: adding ',NbIPtF,'F-points to the',NbIPt,'T-points!')
        zLatLon = np.concatenate( [ zLatLon, zLatLonF ])
        
    return zLatLon





