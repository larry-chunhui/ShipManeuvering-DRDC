from __future__ import division, print_function
import argparse
import numpy as np
from numpy import linspace, zeros, ones, sin, cos, arctan, pi
import os
# Cylinder Re=27000, multiphase flow, 3d
#x streamwise, 
#y perpendicular direction
#z gravity direction 
#input parameters
L=76.44 #Length of wave
alpha_deg=20
theta_deg=45
r=8
R=24
zmax=0.5*L
zmin=-L
xmax=3*L
xmin=-3*L
ymax=L
ymin=-L

print("zmax=",zmax,"zmin=",zmin)
nx=30
nz=100

nx1=1*nx
nx2=2*nx
nx3=3*nx
nx5=5*nx


ny1=nx1
ny2=nx2
ny3=nx3
ny5=nx5

alpha = np.deg2rad(alpha_deg)  # Angle of attack (in radians)
theta = np.deg2rad(theta_deg)

scale = 1            # Scaling factor

# stretch in ring layer
ratio=2
Sn=R-r
q=np.power(ratio,1/(nx2-1))
y=Sn*(1-q)/(1-np.power(q,nx2))
L_theta=2*np.pi*r/(8*nx1)
print("q=",q)
print("y=",y, "L_theta=",L_theta,"dz=",(zmax-zmin)/nz)

# Expansion rates
ExpT = 50          # Expansion rate in transverse direction
ExpD = 10           # Expansion rate in the downstream direction
ExpArc = 5          # Expansion rate along the inlet arc

# ------------------------- END OF MESH PARAMETER REGION --------------------- #


# Create a vector with x-coordinates, camber and thickness


#print(m,p,t)
# Calculate thickness
# The upper expression will give the airfoil a finite thickness at the trailing
# edge, witch might cause trouble. The lower expression is corrected to give
# zero thickness at the trailing edge, but the foil is strictly speaking no
# longer a proper NACA airfoil.
#

print("chunhui")
print(R*cos(theta))
print(np.rad2deg(arctan(1.0)))
print(np.rad2deg(arctan(0.939693/0.34202)))
print(np.rad2deg(arctan(30)))
# Rotate foil to reach specified angle of attack
# Move point of mesh "nose"
# Calculate the location of the vertices on the positive y-axis and put them in a matrix
vertices = zeros((19, 3))

vertices[0, :] = [r,0,zmin]
vertices[1, :] = [R,0,zmin]
vertices[2, :] = [xmax,0,zmin]
vertices[3, :] = [xmax,R*cos(theta),zmin]
vertices[4, :] = [R*cos(theta),R*cos(theta),zmin]
vertices[5, :] = [r*cos(theta),r*cos(theta),zmin]

vertices[6, :] = [xmax,ymax,zmin]
vertices[7, :] = [R*cos(theta),ymax,zmin]
vertices[8, :] = [0,ymax,zmin]
vertices[9, :] = [0,R,zmin]
vertices[10, :] = [0,r,zmin]
vertices[11, :] = [-r,0,zmin]

vertices[12, :] = [-R,0,zmin]
vertices[13, :] = [xmin,0,zmin]
vertices[14, :] = [xmin,R*cos(theta),zmin]

vertices[15, :] = [-R*cos(theta),R*cos(theta),zmin]
vertices[16, :] = [-r*cos(theta),r*cos(theta),zmin]
vertices[17, :] = [xmin,ymax,zmin]
vertices[18, :] = [-R*cos(theta),ymax,zmin]
# Create vertices for other side (negative y-axis)
vertices2 = vertices.copy()
vertices2[:, 2] = zmax
vertices = np.vstack((vertices, vertices2))

vertices3 = zeros((13, 3))
vertices3[0,  :] = [xmax,-R*cos(theta),zmin]
vertices3[1,  :] = [R*cos(theta),-R*cos(theta),zmin]
vertices3[2,  :] = [r*cos(theta),-r*cos(theta),zmin]
vertices3[3,  :] = [xmax,ymin,zmin]
vertices3[4,  :] = [R*cos(theta),ymin,zmin]
vertices3[5,  :] = [0,ymin,zmin]
vertices3[6,  :] = [0,-R,zmin]
vertices3[7,  :] = [0,-r,zmin]

vertices3[8,  :] = [xmin,-R*cos(theta),zmin]
vertices3[9,  :] = [-R*cos(theta),-R*cos(theta),zmin]
vertices3[10,  :] = [-r*cos(theta),-r*cos(theta),zmin]
vertices3[11,  :] = [xmin,ymin,zmin]
vertices3[12,  :] = [-R*cos(theta),ymin,zmin]

vertices4 = vertices3.copy()
vertices4[:, 2] = zmax
vertices3 = np.vstack((vertices3, vertices4))

vertices = np.vstack((vertices, vertices3))
#up
# Edge 0-5 
pts1 = np.array([r*cos(alpha),r*sin(alpha),zmin])
#edge 5-10
pts2 = np.array([r*sin(alpha),r*cos(alpha),zmin])
#edge 16-11
pts3 = np.array([-r*cos(alpha),r*sin(alpha),zmin])
#edge 10-16
pts4 = np.array([-r*sin(alpha),r*cos(alpha),zmin])

# Edge 19-24
pts5 = np.array([r*cos(alpha),r*sin(alpha),zmax])
#edge 24-29
pts6 = np.array([r*sin(alpha),r*cos(alpha),zmax])
#edge 30-35
pts7 = np.array([-r*cos(alpha),r*sin(alpha),zmax])
#edge 35-29
pts8 = np.array([-r*sin(alpha),r*cos(alpha),zmax])

#edge 9-15
pts9 = np.array([-R*sin(alpha),R*cos(alpha),zmin])
#edge 15-12
pts10 = np.array([-R*cos(alpha),R*sin(alpha),zmin])
# Edge 1-4 
pts11 = np.array([R*cos(alpha),R*sin(alpha),zmin])
#edge 4-9
pts12 = np.array([R*sin(alpha),R*cos(alpha),zmin])

#edge 34-28
pts13 = np.array([-R*sin(alpha),R*cos(alpha),zmax])
#edge 31-34
pts14 = np.array([-R*cos(alpha),R*sin(alpha),zmax])
# Edge 20-23
pts15 = np.array([R*cos(alpha),R*sin(alpha),zmax])
#edge 23-28
pts16 = np.array([R*sin(alpha),R*cos(alpha),zmax])



# down
# Edge 0-40 
pts17 = np.array([r*cos(alpha),-r*sin(alpha),zmin])
#edge 40-45
pts18 = np.array([r*sin(alpha),-r*cos(alpha),zmin])
#edge 11-48
pts19 = np.array([-r*cos(alpha),-r*sin(alpha),zmin])
#edge 48-45
pts20 = np.array([-r*sin(alpha),-r*cos(alpha),zmin])

# Edge 19-53
pts21 = np.array([r*cos(alpha),-r*sin(alpha),zmax])
#edge 53-58
pts22 = np.array([r*sin(alpha),-r*cos(alpha),zmax])
#edge 30-61
pts23 = np.array([-r*cos(alpha),-r*sin(alpha),zmax])
#edge 61-58
pts24 = np.array([-r*sin(alpha),-r*cos(alpha),zmax]) 

# 
#edge 47-44
pts25 = np.array([-R*sin(alpha),-R*cos(alpha),zmin])
#edge 12-47
pts26 = np.array([-R*cos(alpha),-R*sin(alpha),zmin])
# Edge 1-39
pts27 = np.array([R*cos(alpha),-R*sin(alpha),zmin])
#edge 39-44
pts28 = np.array([R*sin(alpha),-R*cos(alpha),zmin])

#edge 60-57
pts29 = np.array([-R*sin(alpha),-R*cos(alpha),zmax])
#edge 31-60
pts30 = np.array([-R*cos(alpha),-R*sin(alpha),zmax])
# Edge 20-52
pts31 = np.array([R*cos(alpha),-R*sin(alpha),zmax])
#edge 52-57
pts32 = np.array([R*sin(alpha),-R*cos(alpha),zmax])




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
f.write("//up\n")
#x ridius direction
#y tangential direction
#inner cylinder block number: 0,1, 5,6, 10,11, 15,16
#nx2: number of points r ring, total stretch ratio=2
#4*ny1: number of points theta ring, uniform mesh
f.write("hex (5 4 9 10 24 23 28 29)   (%i %i %i) simpleGrading (2 1 1) \n" % (nx2,ny1,nz))
f.write("hex (0 1 4 5 19 20 23 24)  (%i %i %i) simpleGrading (2 1 1) \n" % (nx2,ny1,nz))
f.write("hex (1 2 3 4 20 21 22 23)  (%i %i %i) simpleGrading (6 1 1) \n" % (nx5,ny1,nz))
f.write("hex (4 3 6 7 23 22 25 26)  (%i %i %i) simpleGrading (6 7 1) \n" % (nx5,ny3,nz))
f.write("hex (9 4 7 8 28 23 26 27) (%i %i %i) simpleGrading (1 7 1) \n" % (nx1,ny3,nz))

f.write("hex (15 16 10 9 34 35 29 28) (%i %i %i) simpleGrading (0.5 1 1) \n" % (nx2,ny1,nz))
f.write("hex (12 11 16 15 31 30 35 34)  (%i %i %i) simpleGrading (0.5 1 1) \n" % (nx2,ny1,nz))
f.write("hex (13 12 15 14 32 31 34 33)  (%i %i %i) simpleGrading (0.15 1 1) \n"% (nx3,ny1,nz))
f.write("hex (14 15 18 17 33 34 37 36) (%i %i %i) simpleGrading (0.15 7 1) \n" % (nx3,ny3,nz))
f.write("hex (15 9 8 18 34 28 27 37) (%i %i %i) simpleGrading (1 7 1) \n"      % (nx1,ny3,nz))

f.write("//dowm \n")
f.write("hex (40 45 44 39 53 58 57 52)   (%i %i %i) simpleGrading (1 2 1) \n" % (nx1,ny2,nz))
f.write("hex (0 40 39 1 19 53 52 20)  (%i %i %i) simpleGrading (1 2 1) \n" % (nx1,ny2,nz))
f.write("hex (1 39 38 2 20 52 51 21)  (%i %i %i) simpleGrading (1 6 1) \n" % (nx1,ny5,nz))
f.write("hex (39 42 41 38 52 55 54 51)  (%i %i %i) simpleGrading (7 6 1) \n" % (nx3,ny5,nz))
f.write("hex (44 43 42 39 57 56 55 52) (%i %i %i) simpleGrading (7 1 1) \n" % (nx3,ny1,nz))

f.write("hex (47 44 45 48 60 57 58 61) (%i %i %i) simpleGrading (1 0.5 1) \n" % (nx1,ny2,nz))
f.write("hex (12 47 48 11 31 60 61 30)  (%i %i %i) simpleGrading (1 0.5 1) \n" % (nx1,ny2,nz))
f.write("hex (13 46 47 12 32 59 60 31)  (%i %i %i) simpleGrading (1 0.15 1) \n"% (nx1,ny3,nz))
f.write("hex (46 49 50 47 59 62 63 60) (%i %i %i) simpleGrading (7 0.15 1) \n" % (nx3,ny3,nz))
f.write("hex (47 50 43 44 60 63 56 57) (%i %i %i) simpleGrading (7 1 1) \n"    % (nx3,ny1,nz))

f.write("); \n")
f.write("\n")

f.write("edges \n")
f.write("( \n")
f.write("//up \n")
f.write("    arc 0 5 (%f %f %f) \n" % tuple(pts1))
f.write("    arc 5 10 (%f %f %f) \n" % tuple(pts2))
f.write("    arc 11 16 (%f %f %f) \n" % tuple(pts3))
f.write("    arc 16 10 (%f %f %f) \n" % tuple(pts4))

f.write("    arc 19 24 (%f %f %f) \n" % tuple(pts5))
f.write("    arc 24 29 (%f %f %f) \n" % tuple(pts6))
f.write("    arc 30 35 (%f %f %f) \n" % tuple(pts7))
f.write("    arc 35 29 (%f %f %f) \n" % tuple(pts8))

f.write("    arc 15 9 (%f %f %f) \n" % tuple(pts9))
f.write("    arc 12 15 (%f %f %f) \n" % tuple(pts10))
f.write("    arc 1 4 (%f %f %f) \n" % tuple(pts11))
f.write("    arc 4 9 (%f %f %f) \n" % tuple(pts12))

f.write("    arc 34 28 (%f %f %f) \n" % tuple(pts13))
f.write("    arc 31 34 (%f %f %f) \n" % tuple(pts14))
f.write("    arc 20 23 (%f %f %f) \n" % tuple(pts15))
f.write("    arc 23 28 (%f %f %f) \n" % tuple(pts16))

f.write(" //down \n")
f.write("    arc 0 40 (%f %f %f) \n"  % tuple(pts17))
f.write("    arc 40 45 (%f %f %f) \n" % tuple(pts18))
f.write("    arc 11 48 (%f %f %f) \n" % tuple(pts19))
f.write("    arc 48 45 (%f %f %f) \n" % tuple(pts20))

f.write("    arc 19 53 (%f %f %f) \n" % tuple(pts21))
f.write("    arc 53 58 (%f %f %f) \n" % tuple(pts22))
f.write("    arc 30 61 (%f %f %f) \n" % tuple(pts23))
f.write("    arc 61 58 (%f %f %f) \n" % tuple(pts24))

f.write("    arc 47 44 (%f %f %f) \n" % tuple(pts25))
f.write("    arc 12 47 (%f %f %f) \n" % tuple(pts26))
f.write("    arc 1 39 (%f %f %f) \n"  % tuple(pts27))
f.write("    arc 39 44 (%f %f %f) \n" % tuple(pts28))

f.write("    arc 60 57 (%f %f %f) \n" % tuple(pts29))
f.write("    arc 31 60 (%f %f %f) \n" % tuple(pts30))
f.write("    arc 20 52 (%f %f %f) \n" % tuple(pts31))
f.write("    arc 52 57 (%f %f %f) \n" % tuple(pts32))

f.write("); \n")
f.write("\n")
f.write("patches \n")
f.write("( \n")

f.write("patch outlet \n")
f.write("        ( \n")
f.write("            (2 3 22 21) \n")
f.write("            (3 6 25 22) \n")
f.write("            (38 2 21 51) \n")
f.write("            (41 38 51 54) \n")
f.write("        ) \n")
f.write("\n")
#rest
f.write("symmetryPlane sym1 \n")
f.write("        ( \n")
f.write("            (7 8 27 26) \n")
f.write("            (6 7 26 25) \n")
f.write("            (8 18 37 27) \n")
f.write("            (18 17 36 37) \n")
f.write("        ) \n")
f.write("\n")

f.write("symmetryPlane sym2 \n")
f.write("        ( \n")
f.write("            (43 42 55 56)  \n")
f.write("            (42 41 54 55)  \n")
f.write("            (50 43 56 63)  \n")
f.write("            (49 50 63 62)  \n")
f.write("        ) \n")
f.write("\n")

f.write("patch inlet \n")
f.write("        ( \n")
f.write("            (14 13 32 33) \n")
f.write("            (17 14 33 36) \n")
f.write("            (13 46 59 32)  \n")
f.write("            (46 49 62 59)  \n")
f.write("        ) \n")
f.write("\n")

f.write("wall cylinder \n")
f.write("        ( \n")
f.write("            (10 5 24 29) \n")
f.write("            (5 0 19 24) \n")
f.write("            (16 10 29 35) \n")
f.write("            (11 16 35 30) \n")

f.write("            (40 45 58 53)  \n")
f.write("             (0 40 53 19)   \n")
f.write("             (45 48 61 58)  \n")
f.write("            (48 11 30 61)  \n")
f.write("        ) \n")


#rest
f.write("patch bottom \n")
f.write("        ( \n")
f.write("           (5 10 9 4) \n \
       (0 5 4 1)  \n \
       (1 4 3 2) \n \
       (4 7 6 3) \n \
       (9 8 7 4) \n \
       (16 15 9 10) \n \
       (11 12 15 16) \n \
       (12 13 14 15) \n \
       (14 17 18 15) \n \
       (15 18 8 9) \n \
       //down  \n \
       (45 40 39 44)  \n \
       (40 0 1 39)  \n \
       (1 2 38 39) \n \
       (39 38 41 42) \n \
       (44 39 42 43) \n \
       (48 45 44 47) \n \
       (11 48 47 12) \n \
       (13 12 47 46) \n \
       (46 47 50 49) \n \
       (47 44 43 50) \n")
f.write("        ) \n")
f.write("\n")


f.write("patch atmosphere \n")
f.write("        ( \n")
f.write("     //up \n \
       (24 23 28 29)\n \
       (20 23 24 19)\n \
       (21 22 23 20)\n \
       (22 25 26 23)\n \
       (26 27 28 23)\n \
       (28 34 35 29)\n \
       (34 31 30 35)\n \
       (33 32 31 34)\n \
       (36 33 34 37)\n \
       (37 34 28 27)\n \
       //down\n \
       (53 58 57 52)\n \
       (19 53 52 20)\n \
       (21 20 52 51)\n \
       (51 52 55 54)\n \
       (52 57 56 55)\n \
       (58 61 60 57)\n \
       (61 30 31 60)\n \
       (31 32 59 60)\n \
       (60 59 62 63)\n \
       (57 60 63 56) \n")
f.write("        ) \n")
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