import numpy as np
from vector_class import Vector
from panel_method_class import panel_velocity2,panel_velocity,panel_velocity3
from matplotlib import pyplot as plt
from plot_functions import plot_savedCp_SurfaceContours
from mesh_class import PanelMesh
from LSQ import lsq


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
mu = []                         ## 先在集群上计算出mu的值再回来计算Cp
with open("mu.txt","r") as file:
    for line in file:
        mu.append(float(line))
print(len(mu))
print(len(shells))
mu = mu[:len(shells)]
for j,m in enumerate(mu):
    mesh.panels[j].mu = m

for i,panel in enumerate(mesh.panels):
            
    panel_neighbours = mesh.give_neighbours3(panel,3)
    panel.sigma = - (panel.n * V_fs)
    panel.Velocity = panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2)
            
    panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2

Cp = [panel.Cp for panel in mesh.panels]
plot_savedCp_SurfaceContours(mesh.panels, Cp, elevation=20, azimuth=45)
