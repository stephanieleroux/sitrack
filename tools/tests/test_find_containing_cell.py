#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

from sys import argv, exit
from os import path, mkdir
import numpy as np
from re import split
from netCDF4 import Dataset
#from math import atan2,pi
import mojito   as mjt
#toDegrees = 180./pi



if __name__ == '__main__':


    if not len(argv) in [3]:
        print('Usage: '+argv[0]+'  <mesh_mask> <lat,lon>')
        exit(0)

    cf_mm   = argv[1]
    clatlon = argv[2]

    vll = split( ',',clatlon)

    rlat, rlon = float(vll[0]), float(vll[1])

    print('\n * Chosen target coordinate point:', rlat, rlon)

    
    # Getting model grid metrics and friends:
    mjt.chck4f( cf_mm )

    # Reading mesh metrics into mesh-mask file:
    with Dataset(cf_mm) as id_mm:
        #maskT = id_mm.variables['tmask'][0,0,:,:]
        xlonT  = id_mm.variables['glamt'][0,:,:]
        xlatT  = id_mm.variables['gphit'][0,:,:]        
        xlonF  = id_mm.variables['glamf'][0,:,:]
        xlatF  = id_mm.variables['gphif'][0,:,:]

    #(nj,ni) = np.shape(maskT)
    #maskT = np.array(maskT, dtype='i1')


    [jN,iN], kVert = mjt.FCC( (rlat,rlon), xlatT, xlonT, xlatF, xlonF, cellType='T', rd_found_km=10., iverbose=2 )


    print('\n * We found:')
    print('   => indices for points at center of cell:', jN,iN )

    print('   => indices of the 4 corner points (vertices) of cell (anti-clockwize starting llc):\n',kVert)


    ik = mjt.Plot1Mesh( (rlat,rlon), xlatF, xlonF, kVert, vnames=['llc','lrc','urc','ulc'],
                        pcoor_cntr=(xlatT[jN,iN], xlonT[jN,iN]),
                        fig_name='mesh.png', pcoor_extra=(-999.,-999.), label_extra=None )

    
