from __future__ import division, print_function
import argparse
import numpy as np
from numpy import linspace, zeros, ones, sin, cos, arctan, pi
import os

#input a1,q,Sn
#output a_n,n
def refine(a_1,q,Sn):
    a_n=(a_1-(1.0-q)*Sn)/q
    n=int(np.log(a_n/a_1)/np.log(q)+1)
    return(a_n,n)
# Square Cylinder Re=200, single phase flow, 2d
#x streamwise, 
#y perpendicular direction
#z gravity direction 
#input parameters
stretch_ratio=1.1
B=1
zmax= 2
zmin=-2
xmax= 9
xmin=-17*B
ymax= 10*B
ymin=-10*B

dx_min=0.1
dy_min=dx_min
nx=int(B/dx_min)
ny=nx

nz=10

scale = 1            # Scaling factor

#upstreams
one=1.0
Lx_up=xmax-B/2
Ly_up=ymax-B/2
Lx_down=-B/2-xmin

print(Lx_up, Lx_down, Ly_up)
dx_up,nx_up=refine(dx_min,stretch_ratio,Lx_up)
ratio_x_up=dx_up/dx_min

dx_down,nx_down=refine(dx_min,stretch_ratio,Lx_down)
ratio_x_down=dx_min/dx_down


dy_up,ny_up=refine(dy_min,stretch_ratio,Ly_up)
ratio_y_up=dy_up/dy_min
ratio_y_down=1/ratio_y_up
ny_down=ny_up

# ------------------------- END OF MESH PARAMETER REGION --------------------- #
# Create a vector with x-coordinates, camber and thickness
#print(m,p,t)
# Calculate thickness
# The upper expression will give the airfoil a finite thickness at the trailing
# edge, witch might cause trouble. The lower expression is corrected to give
# zero thickness at the trailing edge, but the foil is strictly speaking no
# longer a proper NACA airfoil.
#
# Rotate foil to reach specified angle of attack
# Move point of mesh "nose"
# Calculate the location of the vertices on the positive y-axis and put them in a matrix
vertices = zeros((16, 3))

vertices[0, :] = [ B/2, B/2,zmin]
vertices[1, :] = [ B/2,-B/2,zmin]
vertices[2, :] = [-B/2,-B/2,zmin]
vertices[3, :] = [-B/2, B/2,zmin]

vertices[4, :] = [xmax,ymax,zmin]
vertices[5, :] = [xmax,ymin,zmin]
vertices[6, :] = [xmin,ymin,zmin]
vertices[7, :] = [xmin,ymax,zmin]

vertices[8, :] = [-B/2,ymax,zmin]
vertices[9, :] = [ B/2,ymax,zmin]
vertices[10, :] =[xmax, B/2,zmin]
vertices[11, :] =[xmax,-B/2,zmin]

vertices[12, :] = [ B/2,ymin,zmin]
vertices[13, :] = [-B/2,ymin,zmin]
vertices[14, :] = [xmin,-B/2,zmin]
vertices[15, :] = [xmin, B/2,zmin]

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
f.write("|  \\\\    /   O peration     | Version:  3.0.x                                 | \n")
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
#x ridius direction
#y tangential direction
#inner cylinder block number: 0,1, 5,6, 10,11, 15,16
#nx2: number of points r ring, total stretch ratio=2
#4*ny1: number of points theta ring, uniform mesh
#nx_down=10
#nx_up=10
#nx=10
#ny_up=10
#ny_down=10
#nz=1
#ratio_x_down=1
#ratio_x_up=ratio_x_down
#ratio_y_up=ratio_x_down
#ratio_y_down=ratio_x_down
f.write("hex ( 15 3 8 7  31 19 24 23)     (%i %i %i) simpleGrading (%f %f 1) \n" % (nx_down,ny_up,  nz, ratio_x_down,ratio_y_up))
f.write("hex ( 14  2 3 15   30 18 19 31)  (%i %i %i) simpleGrading (%f %f 1) \n" % (nx_down,ny,     nz, ratio_x_down,one))
f.write("hex ( 6 13 2 14   22 29 18 30)   (%i %i %i) simpleGrading (%f %f 1) \n" % (nx_down,ny_down,nz, ratio_x_down,ratio_y_down))
f.write("hex ( 13 12  1 2  29 28 17 18 )  (%i %i %i) simpleGrading (%f %f 1) \n" % (nx,ny_down,  nz,         one,  ratio_y_down))
f.write("hex (  12 5 11 1   28 21 27 17)  (%i %i %i) simpleGrading (%f %f 1) \n" % (nx_up,ny_down,nz,   ratio_x_up,ratio_y_down))
f.write("hex ( 1 11 10  0   17 27 26 16)  (%i %i %i) simpleGrading (%f %f 1) \n" % (nx_up,ny,nz,        ratio_x_up,one))
f.write("hex ( 0 10 4 9    16 26 20 25 )  (%i %i %i) simpleGrading (%f %f 1) \n" % (nx_up,ny_up,nz,     ratio_x_up,ratio_y_up))
f.write("hex ( 3 0 9 8     19 16 25 24)  (%i %i %i) simpleGrading (%f %f 1) \n"%  (nx,ny_up,nz,        one,ratio_y_up))



f.write("); \n")
f.write("\n")

f.write("edges \n")
f.write("( \n")
f.write("); \n")
f.write("\n")

f.write("patches \n")
f.write("( \n")

f.write("patch outlet \n")
f.write("        ( \n")
f.write("            (7  15 31 23) \n")
f.write("            (14 15 30 31) \n")
f.write("            (14  6 22 30) \n")
f.write("        ) \n")
f.write("\n")
#rest
f.write("symmetryPlane sym1 \n")
f.write("        ( \n")
f.write("            (4 9 25 20) \n")
f.write("            (9 8 24 25) \n")
f.write("            (8 7 23 24) \n")
f.write("        ) \n")
f.write("\n")

f.write("symmetryPlane sym2 \n")
f.write("        ( \n")
f.write("            (5  12 28 21)  \n")
f.write("            (12 13 29 28)  \n")
f.write("            (13 6  22 29)  \n")
f.write("        ) \n")
f.write("\n")

f.write("patch inlet \n")
f.write("        ( \n")
f.write("            (4  10 26 20) \n")
f.write("            (10 11 27 26) \n")
f.write("            (11 5  21 27) \n")
f.write("        ) \n")
f.write("\n")

f.write("wall cylinder \n")
f.write("        ( \n")
f.write("            (0 1 17 16) \n")
f.write("            (1 2 18 17) \n")
f.write("            (2 3 19 18) \n")
f.write("            (3 0 16 19) \n")
f.write("        ) \n")

#rest
f.write("\n")
f.write("patch bottom \n")
f.write("        ( \n")
f.write("            (15 3 8 7) \n ")
f.write("            (14 2 3 15)\n ")
f.write("            (6 13 2 14)\n ")
f.write("            (13 12 1 2)\n ")
f.write("            (12 5 11 1)\n ")
f.write("            (1 11 10 0)\n ")
f.write("            (0 10 4 9)\n ")
f.write("            (3 0 9 8) \n")
f.write("        ) \n")
f.write("\n")


f.write("patch atmosphere \n")
f.write("        ( \n")
f.write("           (31 19 24 23)\n ")
f.write("           (30 18 19 31)\n ")
f.write("           (22 29 18 30)\n ")
f.write("           (29 28 17 18)\n ")
f.write("           (28 21 27 17)\n ")
f.write("           (17 27 26 16)\n ")
f.write("           (16 26 20 25)\n ")
f.write("           (19 16 25 24)\n ")
f.write("        ) \n")
f.write("     \n")



f.write("\n")
f.write("); \n")
f.write(" \n")
f.write("mergePatchPairs \n")
f.write("( \n")
f.write("); \n")
f.write(" \n")
f.write("// ************************************************************************* // \n")

# Close file
f.close()