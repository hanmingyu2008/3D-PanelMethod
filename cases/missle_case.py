import numpy as np
from vector_class import Vector
from panel_method_class import Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from plot_functions import plot_savedCp_SurfaceContours
from mesh_class import PanelMesh

nodes = []
shells = []
with open("导弹triangle.obj","r") as file:
    ind = 0
    for line in file:
        if ind == 638:
            pass
        line = line.split()
        if len(line)>0 and line[0] == "v":
            x,y,z = [float(t) for t in line[1:]]
            nodes.append((x,y,z))
        if len(line)>0 and line[0] == "f":
            ix,iy,iz = line[1:]
            ix,_ = [int(t) for t in ix.split("/")]
            iy,_ = [int(t) for t in iy.split("/")]
            iz,_ = [int(t) for t in iz.split("/")]
            if ix*iy*iz == 0:
                print(-1)
            shells.append([ix-1,iy-1,iz-1])
        ind += 1
        print(ind)
print(len(nodes))
print(len(shells))

Cp = []
with open("Cp.txt","r") as file:
    for line in file:
        Cp.append(float(line))
print(len(Cp))

mesh = PanelMesh(nodes, shells)
print(1)

# Surface Contour plot
plot_savedCp_SurfaceContours(mesh.panels, Cp, elevation=20, azimuth=45)