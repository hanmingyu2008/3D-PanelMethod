import numpy as np
from vector_class import Vector
from panel_method_class import panel_velocity2
from matplotlib import pyplot as plt
from plot_functions import plot_Cp_SurfaceContours
from mesh_class import PanelMesh


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


mu = []
with open("mu.txt","r") as file:
    for line in file:
        mu.append(float(line))
print(len(mu))
print(len(shells))

mu = mu[:len(mesh.panels)]
print(len(mu))

for i,panel in enumerate(mesh.panels):
            
    panel.Cp = mu[i]

# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)
print(len(mesh.panels))