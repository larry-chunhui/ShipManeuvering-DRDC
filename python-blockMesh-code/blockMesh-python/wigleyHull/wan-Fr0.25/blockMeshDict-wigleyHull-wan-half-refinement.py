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
# half model of wigley by wan de cheng, drift0
#x streamwise, 
#y perpendicular direction
#z gravity direction 
#input parameters
Re=4.67e6
Fr=0.408
rho=998.8
g=9.8

Lpp=2.5
waterline=0.0


beta_deg=10
U=Fr*np.sqrt(g*Lpp)
nu=U*Lpp/Re

beta = np.deg2rad(beta_deg)
Ux=-U*cos(beta)
Uy=-U*sin(beta)
print("U=%f, Ux=%f, Uy=%f" %(U,Ux,Uy)) 
#ship  region
#x(-1.25,1.25)
#y(-0.125,0.125)
#z(-0.156,0.1)
ship_x_max=1.25
ship_y_max=0.125
ship_z_max=0.1

ship_x_min=-1.25
ship_y_min=-0.125
ship_z_min=-0.156


xmax=1.5*Lpp
xmin=-3.5*Lpp
ymax=2.0*Lpp
ymin=-0*Lpp
zmax=0.0399*Lpp
zmin=-Lpp

Lx=xmax-xmin
Ly=ymax-ymin
Lz=zmax-zmin
#------------------------#
#multigrading
#------------------------#
stretch_ratio=1.2
dz_min=0.0528
dx_min=dz_min
dy_min=dz_min
#z direction
refine_z_min=ship_z_min
refine_z_max=zmax
z_uniform=zmax-refine_z_min
nz_refine=int(z_uniform/dz_min)

#z_up=zmax-refine_z_max
z_down=refine_z_min-zmin
dz_down,nz_down=refine(dz_min,stretch_ratio,z_down)
ratio_z_down=dz_min/dz_down

#x direction
refine_x_min=-0.5*Lpp
refine_x_max= 0.5*Lpp
x_uniform=refine_x_max-refine_x_min
nx_refine=int(x_uniform/dx_min)

x_up=xmax-refine_x_max
x_down=-(xmin-refine_x_min)
dx_up,nx_up=refine(dx_min,stretch_ratio,x_up)
dx_down,nx_down=refine(dx_min,stretch_ratio,x_down)
ratio_x_up=dx_up/dx_min
ratio_x_down=dx_min/dx_down

#y direction
refine_y_max=0.1*Lpp
y_uniform=refine_y_max-ymin
ny_refine=int(refine_y_max/dy_min)

y_up=ymax-refine_y_max
dy_up,ny_up=refine(dy_min,stretch_ratio,y_up)
ratio_y_up=dy_up/dy_min



#------------------------#
#total mesh nx,ny,nz
#------------------------#
nx=nx_refine+nx_up+nx_down
ny=ny_refine+ny_up
nz=nz_refine+nz_down
print("nx=%i,ny=%i,nz=%i,nxnynz=%i" %(nx,ny,nz,nx*ny*nz))
print("zmax=",zmax,"zmin=",zmin)


scale = 1            # Scaling factor



vertices = zeros((4, 3))

vertices[0, :] = [xmin,ymin,zmin]
vertices[1, :] = [xmax,ymin,zmin]
vertices[2, :] = [xmax,ymax,zmin]
vertices[3, :] = [xmin,ymax,zmin]

# Create vertices for other side (negative y-axis)
vertices2 = vertices.copy()
vertices2[:, 2] = zmax
vertices = np.vstack((vertices, vertices2))

f = open("SnappyHexMeshDict", "w")
f.write("refinementBox     \n")
f.write("{ \n")
f.write("   type    searchableBox; \n")
f.write("   min     (%f %f %f); \n"%(xmin,ymin,refine_z_min))
f.write("   max     (%f %f %f); \n"%(xmax,ymax,refine_z_max))
f.write("} \n\n\n")

f.write("refinementWaterline     \n")
f.write("{ \n")
f.write("   type    searchableBox; \n")
f.write("   min     (%f %f %f); \n"%(xmin,ymin,-0.010*Lpp))
f.write("   max     (%f %f %f); \n"%(xmax,ymax, 0.015*Lpp))
f.write("} \n\n\n")

f.write("refinementRegions \n")
f.write("{ \n")

#f.write("  hull \n")
#f.write("  { \n")
#f.write("   mode distance;\n")
#f.write("   levels ((%f 4) (%f 3) (%f 2) (%f 1) );\n"%(dx_min/4,dx_min*0.75,dx_min*1.75,dx_min*2.75))
#f.write("  } \n")

f.write("  refinementBox \n")
f.write("  { \n")
f.write("   mode inside;\n")
f.write("   levels ((1e15 2));\n")
f.write("  } \n")

f.write("  refinementWaterline \n")
f.write("  { \n")
f.write("   mode inside;\n")
f.write("   levels ((1e15 3));\n")
f.write("  } \n")

f.write("} \n\n\n")

f.write("    locationInMesh (%f %f %f);"%(xmax-0.01,ymax-0.01,zmin+0.01))
f.close()

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
f.write("//stretch ratio=%f,dx_min=%f\n"%(stretch_ratio,dx_min))
f.write("// nx1= %i, ny=%i, nz=%i\n" %(nx,ny,nz))

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
f.write("hex (0 1 2 3 4 5 6 7)   (%i %i %i) \n" %(nx,ny,nz)) 
f.write("simpleGrading \n")
f.write("( \n")

f.write("       (     \n")
f.write("             (%f %f %f)\n" %(x_down,   nx_down,  ratio_x_down))
f.write("             (%f %f %f)\n" %(x_uniform,nx_refine,1))
f.write("             (%f %f %f)\n" %(x_up,     nx_up,    ratio_x_up))
f.write("       )     \n")

f.write("       (     \n")
#f.write("             (%f %f %f)\n" %(x_down,   nx_down,  ratio_x_down))
f.write("             (%f %f %f)\n" %(y_uniform,ny_refine,1))
f.write("             (%f %f %f)\n" %(y_up,     ny_up,    ratio_y_up))
f.write("       )     \n")

f.write("       (     \n")
f.write("             (%f %f %f)\n" %(z_down,   nz_down,  ratio_z_down))
f.write("             (%f %f %f)\n" %(z_uniform,nz_refine,1))
#f.write("             (%f %f %f)\n" %(z_up,     nz_up,    ratio_z_up))
f.write("       )     \n")

f.write(") \n")

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

f.write("midPlane \n")
f.write("        { \n")
f.write("            type symmetryPlane;\n")
f.write("            faces \n")
f.write("            ( \n")
f.write("              (0 1 5 4) \n")
f.write("            );\n")
f.write("        } \n")
f.write("\n")

f.write("side \n")
f.write("        { \n")
f.write("            type patch;\n")
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
f.close()