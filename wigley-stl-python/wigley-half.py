# %load http://matplotlib.org/mpl_examples/mplot3d/trisurf3d_demo2.py
#Lpp=1, B=0.1, T=0.0625, ref:
#URANS and DES analysis for a Wigley hull at extreme drift angles
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

# u, v are parameterisation variables
u = (np.linspace(-0.5,0.5, endpoint=True, num=200) * np.ones((200, 1))).flatten()
v = np.repeat(np.linspace(-0.0625, 0.04, endpoint=True, num=200), repeats=200).flatten()

print(type(u),type(v))
# This is the Mobius mapping, taking a (u, v) pair and returning an x, y, z
# triple

L=1
B=0.1
T=0.0625

x = u
z = v
y = 0.5*B*(1-np.power(2*u/L,2))*(1-np.power(v/T,2))*(v<0)+0.5*B*(1-np.power(2*u/L,2))*(v>=0)

print(type(x),type(y))


# Triangulate parameter space to determine the triangles
tri = mtri.Triangulation(u, v)


from stl import mesh

data = np.zeros(len(tri.triangles), dtype=mesh.Mesh.dtype)
mobius_mesh = mesh.Mesh(data, remove_empty_areas=False)
mobius_mesh.x[:] = x[tri.triangles]
mobius_mesh.y[:] = y[tri.triangles]
mobius_mesh.z[:] = z[tri.triangles]
mobius_mesh.save('wigley.stl')


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')

# The triangles in parameter space determine which x, y, z points are
# connected by an edgefrom stl import mesh

data = np.zeros(len(tri.triangles), dtype=mesh.Mesh.dtype)
mobius_mesh = mesh.Mesh(data, remove_empty_areas=False)
mobius_mesh.x[:] = x[tri.triangles]
mobius_mesh.y[:] = y[tri.triangles]
mobius_mesh.z[:] = z[tri.triangles]
mobius_mesh.save('wigley-half-stern.stl')


