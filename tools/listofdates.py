#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 TO DO:

'''

from sys import argv, exit
#from os import path, mkdir, makedirs, environ
#import numpy as np
#from re import split

from mojito import epoch2clock, clock2epoch

if __name__ == '__main__':

    if not len(argv) in [4]:
        print('Usage: '+argv[0]+' <YYYYMMDD0> <YYYYMMDDend> <dt_days>')
        exit(0)

              
    cdateA, cdateB = argv[1], argv[2]    
    dt_sec = int(argv[3])*24*3600



    cdtA = cdateA[0:4]+'-'+cdateA[4:6]+'-'+cdateA[6:8]+'_00:00:00'
    cdtB = cdateB[0:4]+'-'+cdateB[4:6]+'-'+cdateB[6:8]+'_00:00:00'

    #print(' *** From '+cdtA+' to '+cdtB+' !')

    idateA, idateB = clock2epoch(cdtA), clock2epoch(cdtB)
    
            
    nT = int( (idateB+dt_sec - idateA)/dt_sec )

    #print(nT); exit(0)
    
    #idate = idateA
    for jT in range(nT):
        idate = idateA + jT*dt_sec
        cdate = epoch2clock(idate)
        
        cdout = str.replace(cdate, '-', '')
        cdout = str.replace(cdout, '_00:00:00', '')
        
        print(cdout,end=' ')
    
    print('')
        
