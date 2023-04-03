#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:

   *

'''

from sys import argv, exit
#from os import path, mkdir
import numpy as np
#from re import split
from netCDF4 import Dataset
#from math import atan2,pi
#import mojito   as mjt

from cartopy.crs import PlateCarree, NorthPolarStereo


import mojito as mjt

from climporn import dump_2d_field

idebug=0
iplot=1


if __name__ == '__main__':



    if not len(argv) in [2]:
        print('Usage: '+argv[0]+' <mesh_mask>')
        exit(0)

    cf_mm = argv[1]


    # Getting model grid metrics and friends:
    imaskt, xlatT, xlonT, xYt, xXt, xYf, xXf, xResKM = mjt.GetModelGrid( cf_mm )


    # Get extra U,V-point metrics:
    xYv, xXv, xYu, xXu = mjt.GetModelUVGrid( cf_mm )

    # Allocating arrays for model data:
    #(Nj,Ni) = np.shape( imaskt )
    #xUu = np.zeros((Nj,Ni))
    #xVv = np.zeros((Nj,Ni))
    #xIC = np.zeros((Nj,Ni)) ; # Sea-ice concentration



    dump_2d_field( 'latT.nc', xlatT, name='glamT' ) #, unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' )
    dump_2d_field( 'lonT.nc', xlonT, name='gphiT' ) #, unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' )


    zY, zX = mjt.NorthStereoProj( xlatT, xlonT, lam0=-45., phi0=70. )

    dump_2d_field( 'Y.nc', zY, name='Y' ) #, unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' )
    dump_2d_field( 'X.nc', zX, name='X' ) #, unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' )    



    # Cartopy does:

    crs_src = PlateCarree()
    crs_trg = NorthPolarStereo(central_longitude=-45, true_scale_latitude=70)

    XX = crs_trg.transform_points(crs_src, xlonT, xlatT) / 1000.

    #print(np.shape(XX), np.shape(xlonT))

    dump_2d_field( 'YC.nc', XX[:,:,1], name='Y' ) #, unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' )
    dump_2d_field( 'XC.nc', XX[:,:,0], name='X' ) #, unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' )
