import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_method_class import panel_velocity2,panel_velocity,panel_velocity3
from matplotlib import pyplot as plt
from PanelMethod3D.plot_functions import plot_savedCp_SurfaceContours
from PanelMethod3D.mesh_class import PanelMesh
from PanelMethod3D.LSQ import lsq


nodes = []
shells = []

with open("wake_test_fuselage_wake.inp","r") as filee:
    for line in filee:
        if line[0] == " ":
            x,y,z = [float(t) for t in line.split()]
            nodes.append((x,y,z))
        if line[0] == "1" and line[1] != "0":
            labels = [int(t)-1 for t in line[1:].split()]
            shells.append(labels[:4])
        if line[0] == "2":
            labels = [int(t)-1 for t in line[1:].split()]
            shells.append(labels[:3])
mesh = PanelMesh(nodes, shells)

V_fs = Vector((1,0,0))
V_fs_norm = V_fs.norm()
Cp = []                         ## 先在集群上计算出mu的值再回来计算Cp
with open("Cp.txt","r") as file:
    for line in file:
        Cp.append(float(line))
print(len(Cp))
print(len(shells))
mu = Cp[:len(shells)]
L = 0
for j,m in enumerate(Cp):
    mesh.panels[j].Cp = m if m > -0.5 else 0
    if m < -0.5:
        print(mesh.panels[j].area,m)
    L -= m * (mesh.panels[j].n * Vector((0,0,1)))

Cp = [panel.Cp for panel in mesh.panels]
plot_savedCp_SurfaceContours(mesh.panels, Cp, elevation=20, azimuth=45)
print(L)