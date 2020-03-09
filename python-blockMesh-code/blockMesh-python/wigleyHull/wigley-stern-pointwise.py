from __future__ import division, print_function
import argparse
import numpy as np
from numpy import linspace, zeros, ones, sin, cos, arctan, pi
import os
# 
#x streamwise, 
#y perpendicular direction
#z gravity direction 
#input parameters
Re=2.40e6
Fr=0.267
rho=998.8
g=9.8

Lpp=1
beta_deg=10
U=Fr*np.sqrt(g*Lpp)
nu=U*Lpp/Re

beta = np.deg2rad(beta_deg)
Ux=-U*cos(beta)
Uy=-U*sin(beta)
print("U=%f, Ux=%f, Uy=%f" %(U,Ux,Uy)) 



xmax=2*Lpp
xmin=-4*Lpp
ymax=2*Lpp
ymin=0

zmax=0.5*Lpp
zmin=-2*Lpp
z_water=0
z_hull_min=0.0625*Lpp

Lx=xmax-xmin
Ly=ymax-ymin
Lz=zmax-zmin

dx=0.07
dy=dx
dz=dx

nx=int(Lx/dx)
ny=int(Ly/dy)
nz=int(Lz/dz)
print("nx=%i,ny=%i,nz=%i" %(nx,ny,nz))

print("zmax=",zmax,"zmin=",zmin)
scale = 1            # Scaling factor

# Expansion rates
vertices = zeros((4, 3))

vertices[0, :] = [xmax,ymin,zmax] 
vertices[1, :] = [xmax,ymax,zmax]
vertices[2, :] = [xmin,ymax,zmax]
vertices[3, :] = [xmin,ymin,zmax]

# Create vertices for other side (negative y-axis)
vertices2 = vertices.copy()
vertices2[:, 2] = z_water
vertices3 = vertices.copy()
vertices3[:, 2] = z_hull_min
vertices4 = vertices.copy()
vertices4[:, 2] = zmin

vertices = np.vstack((vertices, vertices2))
vertices = np.vstack((vertices, vertices3))
vertices = np.vstack((vertices, vertices4))
# Open file
f = open("pointwise-stern.dat", "w")


# Write file

f.write("vertices \n")
f.write("( \n")
for vertex in vertices:
    f.write("    (%f %f %f)\n" % tuple(vertex))
f.write("); \n")
f.write("\n")

f.close()