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
Re=4.67e6
Fr=0.408
rho=998.8
g=9.8

Lpp=2.5
beta_deg=10
U=Fr*np.sqrt(g*Lpp)
nu=U*Lpp/Re

beta = np.deg2rad(beta_deg)
Ux=-U*cos(beta)
Uy=-U*sin(beta)
print("U=%f, Ux=%f, Uy=%f" %(U,Ux,Uy)) 



xmax=1.5*Lpp
xmin=-3.5*Lpp
ymax=Lpp
ymin=-Lpp
zmax=0.0399*Lpp
zmin=-Lpp

Lx=xmax-xmin
Ly=ymax-ymin
Lz=zmax-zmin

dx=0.1
dy=dx
dz=dx

nx=int(Lx/dx)
ny=int(Ly/dy)
nz=int(Lz/dz)
print("nx=%i,ny=%i,nz=%i" %(nx,ny,nz))


print("zmax=",zmax,"zmin=",zmin)




scale = 1            # Scaling factor

# stretch in ring layer
#ratio=2
#Sn=R-r
#q=np.power(ratio,1/(nx2-1))
#y=Sn*(1-q)/(1-np.power(q,nx2))
#L_theta=2*np.pi*r/(8*nx1)
#print("q=",q)
#print("y=",y, "L_theta=",L_theta,"dz=",(zmax-zmin)/nz)

# Expansion rates


print("chunhui")


vertices = zeros((4, 3))

vertices[0, :] = [xmin,ymin,zmin]
vertices[1, :] = [xmax,ymin,zmin]
vertices[2, :] = [xmax,ymax,zmin]
vertices[3, :] = [xmin,ymax,zmin]

# Create vertices for other side (negative y-axis)
vertices2 = vertices.copy()
vertices2[:, 2] = zmax
vertices = np.vstack((vertices, vertices2))


# Open file
f = open("blockMeshDict", "w")


# Write file
f.write("/*--------------------------------*- C++ -*----------------------------------*\\ \n")
f.write("| =========                 |                                                 | \n")
f.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n")
f.write("|  \\\\    /   O peration     | Version:  1712                                 | \n")
f.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      | \n")
f.write("|    \\\\/     M anipulation  |                                                 | \n")
f.write("\\*---------------------------------------------------------------------------*/ \n")
f.write("FoamFile                                                                        \n")
f.write("{                                                                               \n")
f.write("    version     2.0;                                                            \n")
f.write("    format      ascii;                                                          \n")
f.write("    class       dictionary;                                                     \n")
f.write("    object      blockMeshDict;                                                  \n")
f.write("}                                                                               \n")
f.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n")
f.write("\n")
f.write("convertToMeters %f; \n" % scale)
f.write("\n")
f.write("// nx1= %i, nz=%i\n" %(nx,nz))
f.write("vertices \n")
f.write("( \n")
for vertex in vertices:
    f.write("    (%f %f %f)\n" % tuple(vertex))
f.write("); \n")
f.write("\n")

f.write("blocks \n")
f.write("( \n")
f.write("//up\n")
#x ridius direction
#y tangential direction
#inner cylinder block number: 0,1, 5,6, 10,11, 15,16
#nx2: number of points r ring, total stretch ratio=2
#4*ny1: number of points theta ring, uniform mesh
f.write("hex (0 1 2 3 4 5 6 7)   (%i %i %i) simpleGrading (1 1 1) \n" % (nx,ny,nz))
f.write("); \n")
f.write("\n")

f.write("edges \n")
f.write("( \n")
f.write("); \n")

f.write("\n")
f.write("boundary \n")
f.write("( \n")

f.write("atmosphere \n")
f.write("        { \n")
f.write("            type patch;\n")
f.write("            faces \n")
f.write("            ( \n")
f.write("              (4 5 6 7) \n")
f.write("            );\n")
f.write("        } \n")
f.write("\n")

f.write("inlet \n")
f.write("        { \n")
f.write("            type patch;\n")
f.write("            faces \n")
f.write("            ( \n")
f.write("              (1 2 6 5) \n")
f.write("            );\n")
f.write("        } \n")
f.write("\n")

f.write("outlet \n")
f.write("        { \n")
f.write("            type patch;\n")
f.write("            faces \n")
f.write("            ( \n")
f.write("              (0 4 7 3) \n")
f.write("            );\n")
f.write("        } \n")
f.write("\n")

f.write("bottom \n")
f.write("        { \n")
f.write("            type symmetryPlane;\n")
f.write("            faces \n")
f.write("            ( \n")
f.write("              (0 3 2 1) \n")
f.write("            );\n")
f.write("        } \n")
f.write("\n")

f.write("side \n")
f.write("        { \n")
f.write("            type symmetryPlane;\n")
f.write("            faces \n")
f.write("            ( \n")
f.write("              (0 1 5 4) \n")
f.write("            );\n")
f.write("        } \n")
f.write("\n")

f.write("midPlane \n")
f.write("        { \n")
f.write("            type symmetryPlane;\n")
f.write("            faces \n")
f.write("            ( \n")
f.write("              (3 7 6 2) \n")
f.write("            );\n")
f.write("        } \n")
f.write("\n")

f.write(" );\n")
f.write(" \n")
f.write("mergePatchPairs \n")
f.write("( \n")
f.write("); \n")
f.write(" \n")
f.write("// ************************************************************************* // \n")

# Close file
f.close()