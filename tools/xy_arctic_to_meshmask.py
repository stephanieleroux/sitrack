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
from math import copysign

from netCDF4 import Dataset

#import mojito   as mjt
import sitrack  as sit


#import random
#from random import random, choices


rmasked = -999.


# Arctic projection parameters:
rlon0 = -45.
rlat0 =  70.

idebug=0
iplot=1


# intersection between line(p1, p2) and line(p3, p4)
def intersect( ps1, ps2 ):
    '''
        Returns the cartesian coordinates of the intersection of 2 segments
        Returns None if the 2 segmenst do not intersect
    '''
    (y1,x1, y2,x2) = ps1
    (y3,x3, y4,x4) = ps2
    #
    denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
    #
    if denom == 0: # parallel
        return None
    #
    ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
    #
    if ua < 0 or ua > 1: # out of range
        return None
    #
    ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
    #
    if ub < 0 or ub > 1: # out of range
        return None
    #
    x = x1 + ua * (x2-x1)
    y = y1 + ua * (y2-y1)
    #
    return (y,x)



def __argument_parsing__():
    '''
    ARGUMENT PARSING / USAGE
    '''
    import argparse as ap
    #
    parser = ap.ArgumentParser(description='convert geographical coordinates to mesh / metrics info')
    #
    # Required:
    rqrdNam = parser.add_argument_group('required arguments')
    rqrdNam.add_argument('-i', '--fin' , required=True,        help='input file containg longitude and latitude arrays')
    #
    # Optional:
    parser.add_argument('-x', '--nXpc'  , default='x',          help='name of projection x-coordinate (default="x")')
    parser.add_argument('-y', '--nYpc'  , default='y',          help='name of projection y-coordinate (default="y")')
    parser.add_argument('-f', '--field' , default='siconc',     help='name of 2D field to get land-sea mask from")')
    parser.add_argument('-o', '--fout'  , default='coordinates_meshmask.nc',    help='output file (default="coordinates_meshmask.nc")')
    #
    args = parser.parse_args()
    #
    return args.fin, args.nXpc, args.nYpc, args.field, args.fout






if __name__ == '__main__':

    print('')

    cf_in, cv_x, cv_y, cv_field, cf_out = __argument_parsing__()

    print(' *** Input file: '+cf_in)
    print('                => name for longitude and latitude: "'+cv_x+'", "'+cv_y+'"\n')

    sit.chck4f(cf_in)


    with Dataset(cf_in) as id_in:
        shpX, shpY = id_in.variables[cv_x].shape, id_in.variables[cv_y].shape
    
    l_2d_coordinates = ( len(shpX)==2 and len(shpY)==2 )
     
    if l_2d_coordinates:
        if shpX != shpY :
            print('ERROR: lon and lat variables are 2D and have different shapes!'); exit(0)
            print(' *** Coordinates are 2D, irregular grid!')
        #
        (Ny,Nx) = shpX
        print('FIX ME: irregular 2D X,Y grid!'); exit(0)
        #
    elif ( len(shpX)==1 and len(shpY)==1 ):
        print(' *** X,Y coordinates are 1D! => regular cartesian X,Y grid!')
        (Nx,), (Ny,) = shpX, shpY
        #
    else:
        print('ERROR: could not figure out the shape of coordinates in inpute file...'); exit(0)        
    #
    del shpX, shpY
    print('       ==> domain shape: Ny, Nx =', Ny,Nx,'\n')


    xfield = np.zeros((Ny,Nx))
    
    with Dataset(cf_in) as id_in:
        rfill_val   = id_in.variables[cv_field]._FillValue
        print(' *** Fill/Missing value for '+cv_field+' =', rfill_val)
        
        xfield[:,:] = id_in.variables[cv_field][0,:,:]
    
        if not l_2d_coordinates:
            vy = id_in.variables[cv_y][:]
            vx = id_in.variables[cv_x][:]

            # Making the 1D X,Y coordinates 2D and in km:
            xYkm_t, xXkm_t = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
            for i in range(Nx): xYkm_t[:,i] = 1.e-3 * vy[:]
            for j in range(Ny): xXkm_t[j,:] = 1.e-3 * vx[:]
    
    del vy, vx
        
    #print(xXkm_t); print('\n')
    #print(xYkm_t); print('\n')
    #exit(0)



    ### closing `cf_in`...


    print(' *** Creating land-sea mask based on masked values of field "'+cv_field+'"...')
    imask = np.ones((Ny,Nx),dtype='int')
    #imask[np.where(xfield < rfill_val - copysign(0.001,rfill_val))] = 0
    imask[np.where(xfield[:,:] == rfill_val)] = 0
    print(imask[::100,::100],'\n')

    
    print(' *** Allocating arrays...')
    # Geographic coordinates:
    xlat_t,xlon_t = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    xlat_u,xlon_u = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    xlat_v,xlon_v = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    xlat_f,xlon_f = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    # Cartesian coordinates:
    xYkm_u,xXkm_u = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    xYkm_v,xXkm_v = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    xYkm_f,xXkm_f = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    print('       .... done!\n')
    
    #xYkm_t, xXkm_t = sit.ConvertGeo2CartesianNPSkm( xlat_t, xlon_t,  lat0=rlat0, lon0=rlon0 )


    
    # Constructing U-points:
    xYkm_u[:,0:Nx-1] = 0.5 * ( xYkm_t[:,0:Nx-1] + xYkm_t[:,1:Nx] )
    xXkm_u[:,0:Nx-1] = 0.5 * ( xXkm_t[:,0:Nx-1] + xXkm_t[:,1:Nx] )

    # Constructing V-points:
    xYkm_v[0:Ny-1,:] = 0.5 * ( xYkm_t[0:Ny-1,:] + xYkm_t[1:Ny,:] )
    xXkm_v[0:Ny-1,:] = 0.5 * ( xXkm_t[0:Ny-1,:] + xXkm_t[1:Ny,:] )



    print(' *** Building F-points:')
    for j in range(Ny-1):
        for i in range(Nx-1):

            zs1u = ( xYkm_u[j,i],xXkm_u[j,i] , xYkm_u[j+1,i],xXkm_u[j+1,i] ) ; # S-N segment that connects 2 U-point
            zs1v = ( xYkm_v[j,i],xXkm_v[j,i] , xYkm_v[j,i+1],xXkm_v[j,i+1] ) ; # W-E segment that connects 2 V-point

            zint = intersect( zs1u, zs1v )

            if not zint:
                print('ERROR: problem when building F-points...'); exit(0); # that should really not happen!
            else:
                (xYkm_f[j,i],xXkm_f[j,i]) = zint

        ###
        if j%100 == 0: print('        * j =',j,'/',Ny-1)
    ###
    print('')

    print(' *** Converting back to Geographic coordinates, with lat0='+str(rlat0)+', lon0='+str(rlon0)+' !')
    xlat_t, xlon_t = sit.ConvertCartesianNPSkm2Geo( xYkm_t, xXkm_t,  lat0=rlat0, lon0=rlon0 )
    xlat_u, xlon_u = sit.ConvertCartesianNPSkm2Geo( xYkm_u, xXkm_u,  lat0=rlat0, lon0=rlon0 )
    xlat_v, xlon_v = sit.ConvertCartesianNPSkm2Geo( xYkm_v, xXkm_v,  lat0=rlat0, lon0=rlon0 )
    xlat_f, xlon_f = sit.ConvertCartesianNPSkm2Geo( xYkm_f, xXkm_f,  lat0=rlat0, lon0=rlon0 )


    print(' *** Computing e1t and e2t...')
    xe2t, xe1t = np.zeros((Ny,Nx), dtype=np.double)+rmasked, np.zeros((Ny,Nx), dtype=np.double)+rmasked
    xe2t[1:Ny,:] = 1000. * (xYkm_u[1:Ny,:] - xYkm_u[0:Ny-1,:]) # in m 
    xe1t[:,1:Nx] = 1000. * (xXkm_u[:,1:Nx] - xXkm_u[:,0:Nx-1]) # in m 


    # Writing output file:
    id_out = Dataset(cf_out, 'w', format='NETCDF4')

    id_out.createDimension('y', Ny)
    id_out.createDimension('x', Nx)


    id_mskt  = id_out.createVariable( 'tmask' ,'f4',('y','x',), zlib=True, complevel=7 )
    id_mskt[:,:] = imask[:,:].astype(np.single)
    
    id_lat_t  = id_out.createVariable( 'gphit' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_t[:,:] = xlat_t[:,:] ; id_lat_t.units = 'degrees_north'
    id_lon_t  = id_out.createVariable( 'glamt' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_t[:,:] = xlon_t[:,:] ; id_lon_t.units = 'degrees_east'
    #
    id_lat_u  = id_out.createVariable( 'gphiu' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_u[:,:] = xlat_u[:,:] ; id_lat_u.units = 'degrees_north'
    id_lon_u  = id_out.createVariable( 'glamu' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_u[:,:] = xlon_u[:,:] ; id_lon_u.units = 'degrees_east'    
    #
    id_lat_v  = id_out.createVariable( 'gphiv' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_v[:,:] = xlat_v[:,:] ; id_lat_v.units = 'degrees_north'
    id_lon_v  = id_out.createVariable( 'glamv' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_v[:,:] = xlon_v[:,:] ; id_lon_v.units = 'degrees_east'
    #
    id_lat_f  = id_out.createVariable( 'gphif' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_f[:,:] = xlat_f[:,:] ; id_lat_f.units = 'degrees_north'
    id_lon_f  = id_out.createVariable( 'glamf' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_f[:,:] = xlon_f[:,:] ; id_lon_f.units = 'degrees_east'
    #
    #
    id_Ykm_t  = id_out.createVariable( 'ypos_t' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Ykm_t[:,:] = xYkm_t[:,:] ; id_Ykm_t.units = 'km'
    id_Xkm_t  = id_out.createVariable( 'xpos_t' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Xkm_t[:,:] = xXkm_t[:,:] ; id_Xkm_t.units = 'km'
    #
    id_Ykm_u  = id_out.createVariable( 'ypos_u' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Ykm_u[:,:] = xYkm_u[:,:] ; id_Ykm_u.units = 'km'
    id_Xkm_u  = id_out.createVariable( 'xpos_u' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Xkm_u[:,:] = xXkm_u[:,:] ; id_Xkm_u.units = 'km'    
    #
    id_Ykm_v  = id_out.createVariable( 'ypos_v' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Ykm_v[:,:] = xYkm_v[:,:] ; id_Ykm_v.units = 'km'
    id_Xkm_v  = id_out.createVariable( 'xpos_v' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Xkm_v[:,:] = xXkm_v[:,:] ; id_Xkm_v.units = 'km'
    #
    id_Ykm_f  = id_out.createVariable( 'ypos_f' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Ykm_f[:,:] = xYkm_f[:,:] ; id_Ykm_f.units = 'km'
    id_Xkm_f  = id_out.createVariable( 'xpos_f' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_Xkm_f[:,:] = xXkm_f[:,:] ; id_Xkm_f.units = 'km'
    #
    #
    id_e2t  = id_out.createVariable( 'e2t' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_e2t[:,:] = xe2t[:,:] ; id_e2t.units = 'm'
    id_e1t  = id_out.createVariable( 'e1t' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_e1t[:,:] = xe1t[:,:] ; id_e1t.units = 'm'
    #
    #
    id_out.about = '`lat_lon_to_mesh.py` of `sitrack`, based on file "'+cf_in+'"'
    id_out.close()

    print('\n *** '+cf_out+' written!\n')




    if iplot>0:

        import climporn as cp
        
        cfig = 'map_mesh'
        isubsamp = 10
        DPIsvg = 100
        rLat0 = 75.
        
        ii = cp.PlotGridGlobe( xlon_f[:,:], xlat_f[:,:],
                               chemi='N', lon0=-35., lat0=rLat0, cfig_name=cfig+'_NH_35W_f_OUT_ortho_WHITE.svg',
                               nsubsamp=isubsamp, rdpi=DPIsvg, ldark=False )

        
    



    
