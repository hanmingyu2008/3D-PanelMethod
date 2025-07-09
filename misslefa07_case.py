import numpy as np
from vector_class import Vector
from panel_method_class import Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from plot_functions import plot_Cp_SurfaceContours
from mesh_class import PanelMesh

nodes = []
shells = []
with open("fa07morepoint.obj","r") as file:
    for line in file:
        if line[0] == "#":
            pass
        if line[0] == "v":
            x,y,z = [float(t) for t in line[1:].split()]
            nodes.append((x,y,z))
        if line[0] == "f":
            ix,iy,iz = [int(t) for t in line[1:].split()]
            if ix*iy*iz == 0:
                print(-1)
            shells.append([ix-1,iy-1,iz-1])
print(len(nodes))
print(len(shells))

mesh = PanelMesh(nodes, shells)
print(1)

# mesh = PanelMesh.generate_from_stl_file("Sphere_0080")


V_fs = Vector((1, 0, 0))

panel_method = Steady_Wakeless_PanelMethod(V_fs)
print(2)
panel_method.solve(mesh)
print(3)

# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)