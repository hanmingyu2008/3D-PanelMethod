import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_method_class import panel_velocity2,panel_velocity,panel_velocity3,Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from PanelMethod3D.plot_functions import plot_savedCp_SurfaceContours,plot_Cp_SurfaceContours
from PanelMethod3D.mesh_class import PanelMesh
from PanelMethod3D.LSQ import lsq
from fuckingdatacleaning import find_positive_faces

nodes = []
shells = []
temp = {}
with open("hex.obj","r") as filee:
    for line in filee:
        lst = line.split()
        if len(lst) == 0:
            continue
        if lst[0] == "v":
            x,y,z = [float(t) for t in lst[1:]]
            nodes.append((x,y,z))
        if lst[0] == "f":
            a,b,c = [int(t.split("/")[0])-1 for t in lst[1:]]
            if a!=b and b!=c and c!=a:
                shells.append([a,b,c])

print(len(shells))
lst = find_positive_faces(nodes,shells)
print(lst)
print(len(lst))
mesh = PanelMesh(nodes,shells)
# method = Steady_Wakeless_PanelMethod(Vector((0,0,1)))
# method.solve(mesh)
mesh.panels[lst[0]].Cp = 1
plot_Cp_SurfaceContours([mesh.panels[i] for i in lst])