# This adjustment is a 3d Conformal Adjustment that will be done using a combined Functional Model
# Approximating for 7 parameters: Scale, Rotations, Translations

from data import *
from sympy import *
import pandas as pd
import numpy as np
import math as m
from numpy.linalg import inv

np.set_printoptions(precision=8,threshold=1000, edgeitems=4,linewidth = 800, suppress=True)

# UNKNOWNS: w  p  k  Tx  Ty Tz  s


# Initial Approximations 
# Scale - Distance between points 1 and 2 in both Coordinate Systems:
scale = (sqrt((C2x - C1x)**2 + (C2y - C1y)**2 + (C2z - C1z)**2))/(sqrt((c2x - c1x)**2 + (c2y - c1y)**2 + (c2z - c1z)**2))
scale = scale.subs(variables)

# Rotation Parameters
# First Determine which triangle in XYZ coordinate system has strongest geometery using altitude:
def altitude(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    a = (x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2
    b = (x3 - x2)**2 + (y3 - y2)**2 + (z3 - z2)**2
    c = (x1 - x3)**2 + (y1 - y3)**2 + (z1 - z3)**2
    sides = np.array([a.subs(variables),b.subs(variables),c.subs(variables)])
    sides = -np.sort(-sides)
    h = sides[1] - ((sides[0] + sides[1] - sides[2])/(2*sqrt(sides[0])))**2
    return sqrt(h) 

# Find Triangle with largest Altitude
tri_123 = altitude(C1x,C1y,C1z,C2x,C2y,C2z,C3x,C3y,C3z)
tri_124 = altitude(C1x,C1y,C1z,C2x,C2y,C2z,C4x,C4y,C4z)
tri_143 = altitude(C1x,C1y,C1z,C4x,C4y,C4z,C3x,C3y,C3z)
tri_243 = altitude(C2x,C2y,C2z,C4x,C4y,C4z,C3x,C3y,C3z)

# Find Normal Vectors - Using Points 2,3 and 4
def normal(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    a = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1)
    b = (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1)
    c = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    a = a.subs(variables)
    b = b.subs(variables)
    c = c.subs(variables)
    return np.array([a,b,c])

norm_XYZ = normal(C2x,C2y,C2z,C3x,C3y,C3z,C4x,C4y,C4z)
norm_xyz = normal(c2x,c2y,c2z,c3x,c3y,c3z,c4x,c4y,c4z)


# Find the Tilt and Azimuth of each System
tilt_XYZ = atan(norm_XYZ[2]/sqrt(norm_XYZ[0]**2 + norm_XYZ[1]**2))
tilt_xyz = atan(norm_xyz[2]/sqrt(norm_xyz[0]**2 + norm_xyz[1]**2))
Az_XYZ = atan(norm_XYZ[0]/norm_XYZ[1])
Az_xyz = atan(norm_xyz[0]/norm_xyz[1])

# Determine the Rotation matrix using swing value of 0
def rotation(a,t,s):
    a = m.radians(a) # Azimuth 
    t = m.radians(t) # Tilt
    m11 = - cos(a)*cos(s) - sin(a)*cos(t)*sin(s)
    m12 = sin(a)*cos(s) - cos(a)*cos(t)*sin(s)
    m13 = - sin(t)*sin(s)
    m21 = cos(a)*sin(s) - sin(a)*cos(t)*cos(s)
    m22 = - sin(a)*sin(s) - cos(a)*cos(t)*cos(s)
    m23 = - sin(t)*cos(s)
    m31 = - sin(a)*sin(a)
    m32 = - cos(a)*sin(t)
    m33 = cos(t)
    return np.array([[m11, m12, m13],[m21, m22, m23],[m31, m32, m33]],dtype=float)
    
rot_XYZ = rotation(Az_XYZ,tilt_XYZ,0)
rot_xyz = rotation(Az_xyz,tilt_xyz,0)

# Determine the Swing - Difference of azimuths for common line
# This is done using 2 points in each coordinate system 
# These are the same points as were found to build the strongest triangle (234)
def Az_control(mat,x,y,z,x1,y1,z1):
    xp1 = mat[0,0]*x + mat[0,1]*y + mat[0,2]*z
    yp1 = mat[1,0]*x + mat[1,1]*y + mat[1,2]*z
    xp2 = mat[0,0]*x1 + mat[0,1]*y1 + mat[0,2]*z1
    yp2 = mat[1,0]*x1 + mat[1,1]*y1 + mat[1,2]*z1
    xp1 = xp1.subs(variables)
    xp2 = xp2.subs(variables)
    yp1 = yp1.subs(variables)
    yp2 = yp2.subs(variables)
    return atan(m.radians(xp2 - xp1)/m.radians(yp2 - yp1)) 

Ac_XYZ = Az_control(rot_XYZ,C2x,C2y,C2z,C3x,C3y,C3z)
Ac_xyz = Az_control(rot_xyz,c2x,c2y,c2z,c3x,c3y,c3z)
swing = Ac_XYZ - Ac_xyz

# Overall Rotation Matrix
M1 = rotation(Az_xyz,tilt_xyz,swing)
M11 = M1
M12 = rot_XYZ
M = M11.T@M12

# Approximations for Omega, Phi, Kappa
omega = atan(-M[2,1]/M[2,2])
phi = asin(M[2,0])
kappa = atan(-M[1,0]/M[0,0])

# Approximations of Translations
T_x = (C1x - scale*(c1x*M[0,0] + c1y*M[1,0] + c1z*M[2,0])).subs(variables)
T_y = (C1y - scale*(c1x*M[0,1] + c1y*M[1,1] + c1z*M[2,1])).subs(variables)
T_z = (C1z - scale*(c1x*M[0,2] + c1y*M[1,2] + c1z*M[2,2])).subs(variables)


# First Design Matrix
# Dimensions of 12 x 7
A = np.zeros([12,7],dtype = float)

# Rotations
w,pi,k = symbols('w pi k')
m11 = cos(w)*cos(k)
m12 = sin(w)*sin(pi)*cos(k) + cos(w)*sin(k)
m13 = - cos(w)*sin(pi)*cos(k) + sin(w)*sin(k)
m21 = - cos(pi)*sin(k)
m22 = - sin(w)*sin(pi)*sin(k) + cos(w)*cos(k)
m23 = cos(w)*sin(pi)*sin(k) + sin(w)*cos(k)
m31 = sin(pi)
m32 = - sin(w)*cos(pi)  
m33 = cos(w)*cos(pi)


x,y,z,Tx,Ty,Tz,s = symbols('x y z Tx Ty Tz s')
variables_a = [(w,omega),(pi,phi),(k,kappa),(Tx,T_x),(Ty,T_y),(Tz,T_z),(s,scale)]
var_a = np.array([w,pi,k,Tx,Ty,Tz,s])
# Functional Model - First Design
X = s*(x*m11 + y*m21 + z*m31) + Tx
Y = s*(x*m12 + y*m22 + z*m32) + Ty
Z = s*(x*m13 + y*m23 + z*m33) + Tz

X_a = (diff(X,var_a[:])).subs(variables_a)
Y_a = (diff(Y,var_a[:])).subs(variables_a)
Z_a = (diff(Z,var_a[:])).subs(variables_a)

counter = 0
for i in range(0,12,3):
    x_now = X_a.subs([(x,control_a[i]),(y,control_a[i+1]),(z,control_a[i+2])])
    y_now = Y_a.subs([(x,control_a[i]),(y,control_a[i+1]),(z,control_a[i+2])])
    z_now =Z_a.subs([(x,control_a[i]),(y,control_a[i+1]),(z,control_a[i+2])])
    x_now = x_now.subs(variables)
    y_now = y_now.subs(variables)
    z_now = z_now.subs(variables)
    A[counter,:] = x_now
    A[counter+1,:] = y_now
    A[counter+2,:] = z_now
    counter +=3


# Second Design Matrix
# Math Model Second Design Matrix
Xb, Yb, Zb = symbols('Xb Yb Zb')
Bx = s*(x*m11 + y*m21 + z*m31) - Xb + Tx
By = s*(x*m12 + y*m22 + z*m32) - Yb + Ty
Bz = s*(x*m13 + y*m23 + z*m33) - Zb + Tz

counter = 0
other_counter = 0
B = np.zeros([12,24],dtype=float)
for i in range(0,12,3):
    variables_b = [(x,control_a[counter]),(y,control_a[counter+1]),(z,control_a[counter+2])]
    var_b = np.array([x,y,z,Xb,Yb,Zb])

    X_b = (diff(Bx,var_b[:])).subs(variables_b)
    Y_b = (diff(By,var_b[:])).subs(variables_b)
    Z_b = (diff(Bz,var_b[:])).subs(variables_b)
    x_now = X_b.subs(variables_a)
    y_now = Y_b.subs(variables_a)
    z_now = Z_b.subs(variables_a)
    B[counter,other_counter:other_counter+6] = x_now
    B[counter+1,other_counter:other_counter+6] = y_now
    B[counter+2,other_counter:other_counter+6] = z_now
    counter +=3
    other_counter +=1

# Adjustment
M = B@cl@B.T

# Misclosure w is l - fx
# We can use the coordinates in the arbitrary system as approximates
#of observations in the Control system
#W = l-fxo = XYZ coords - xyz coords 
l = np.array([C1x, C1y, C1z, C2x, C2y, C2z, C3x, C3y, C3z, C4x, C4y, C4z])
Fx = np.array([c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z, c4x, c4y, c4z])
w = l - Fx
for i in range(w.size):
    w[i] = w[i].subs(variables)
w = np.array(w,dtype=float)

# Delta X
dx = inv(A.T@inv(M)@A)@(A.T@inv(M)@w)
print(dx)











