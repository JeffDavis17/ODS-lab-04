import pandas as pd
import numpy as np
from sympy import *


d = pd.read_csv(r'C:/Users/jcdav/Documents/GitHub/ODS-lab-04/code/data.csv')
d = np.array(d)

# Variables
# Common Points
com_pt1 = np.array(d[11:15,1:4],dtype=float) # Common points from Coordinate system 1
com_pt2 = np.array(d[11:15,4:7],dtype=float) # Common points from Coordinate system 2

# Covariances common points
cov_com_pt1 = np.array(d[19:31,1:4],dtype=float)
cov_com_pt2 = np.array(d[19:31,4:7],dtype=float)

# Covariance matrix of observables
cl = np.zeros([24,24],dtype=float)
cl[0:3,0:3] = cov_com_pt2[0:3,0:3]
cl[3:6,3:6] = cov_com_pt1[0:3,0:3]
cl[6:9,6:9] = cov_com_pt2[3:6,0:3]
cl[9:12,9:12] = cov_com_pt1[3:6,0:3]
cl[12:15,12:15] = cov_com_pt2[6:9,0:3]
cl[15:18,15:18] = cov_com_pt1[6:9,0:3]
cl[18:21,18:21] = cov_com_pt2[9:12,0:3]
cl[21:24,21:24] = cov_com_pt1[9:12,0:3]

# Coordinates to be Transformed
p = np.array(d[34:37,1:4],dtype=float)

# Covariances pts to be transformed
cov_p = np.array(d[41:50,1:4],dtype=float)

# Sympy Variables
C1x, C1y, C1z, C2x, C2y, C2z, C3x, C3y, C3z, C4x, C4y, C4z = symbols('C1x, C1y, C1z, C2x, C2y, C2z, C3x, C3y, C3z, C4x, C4y, C4z') 
c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z, c4x, c4y, c4z = symbols('c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z, c4x, c4y, c4z') 
p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z = symbols('p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z')

control_a = [c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z, c4x, c4y, c4z]

var = [(C1x,com_pt1[0,0]),
    (C1y,com_pt1[0,1]),
    (C1z,com_pt1[0,2]),
    (C2x,com_pt1[1,0]),
    (C2y,com_pt1[1,1]),
    (C2z,com_pt1[1,2]),
    (C3x,com_pt1[2,0]),
    (C3y,com_pt1[2,1]),
    (C3z,com_pt1[2,2]),
    (C4x,com_pt1[3,0]),
    (C4y,com_pt1[3,1]),
    (C4z,com_pt1[3,2]),
    (c1x,com_pt2[0,0]),
    (c1y,com_pt2[0,1]),
    (c1z,com_pt2[0,2]),
    (c2x,com_pt2[1,0]),
    (c2y,com_pt2[1,1]),
    (c2z,com_pt2[1,2]),
    (c3x,com_pt2[2,0]),
    (c3y,com_pt2[2,1]),
    (c3z,com_pt2[2,2]),
    (c4x,com_pt2[3,0]),
    (c4y,com_pt2[3,1]),
    (c4z,com_pt2[3,2]),
    (p1x,p[0,0]),
    (p1y,p[0,1]),
    (p1z,p[0,2]),
    (p2x,p[1,0]),
    (p2y,p[1,1]),
    (p2z,p[1,2]),
    (p3x,p[2,0]),
    (p3y,p[2,1]),
    (p3z,p[2,2])]

variables = var








