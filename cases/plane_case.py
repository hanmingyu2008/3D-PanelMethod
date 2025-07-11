import numpy as np
from vector_class import Vector
from panel_method_class import Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from plot_functions import plot_savedCp_SurfaceContours
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
print(len(nodes))
print(len(shells))
print(len(nodes))
print(len(shells))

Cp = []
with open("Cp.txt","r") as file:
    for line in file:
        Cp.append(float(line))
print(len(Cp))

mesh = PanelMesh(nodes, shells)
Cp = Cp[:len(mesh.panels)]
print(len(Cp))
print(1)

# Surface Contour plot
plot_savedCp_SurfaceContours(mesh.panels, Cp, elevation=20, azimuth=45)
print(len(mesh.panels))