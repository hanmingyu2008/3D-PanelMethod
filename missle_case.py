import numpy as np
from vector_class import Vector
from panel_method_class import panel_velocity2,panel_velocity,panel_velocity3
from matplotlib import pyplot as plt
from plot_functions import plot_savedCp_SurfaceContours,plot_Cp_SurfaceContours
from mesh_class import PanelMesh
from LSQ import lsq

nodes = []
shells = []
with open("model_normalized.obj","r") as filee:
    for line in filee:
        lst = line.split()
        if len(lst) == 0:
            continue
        if lst[0] == "v":
            x,y,z = [float(t) for t in lst[1:]]
            nodes.append((x,y,z))
        if lst[0] == "f":
            a,b,c = [int(t.split("/")[0])-1 for t in lst[1:]]
            shells.append([a,b,c])
print(1)
mesh = PanelMesh(nodes,shells)
print(len(mesh.panels))
print(mesh.node_num)
mesh.panels[0].Cp = 1
plot_Cp_SurfaceContours(mesh.panels)