#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################

from sys import argv, exit
#from os import path
import numpy as np

#from scipy.spatial import Delaunay

import mojito   as mjt


# Quad as [ [x0,y0], [x1,y1], [x2,y2], [x3,y3] ]
quad1 = np.array([ [0.,0.], [3.,0.], [4.,4.], [1.,3.5] ])


test_points = np.array( [
    [2.,2.],
    [6.,6.],
    [-1.,2.],    
    [3.1,3.6],
])

(ntP,_) = np.shape(test_points)

for jtP in range(ntP):
    print('Inside',test_points[jtP],'? =>',mjt.IsInsideQuadrangle( test_points[jtP,0], test_points[jtP,1], quad1 ))
