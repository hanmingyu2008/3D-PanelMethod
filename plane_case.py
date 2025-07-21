import numpy as np
from vector_class import Vector
from panel_method_class import panel_velocity2,panel_velocity,panel_velocity3
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
mesh = PanelMesh(nodes, shells)
'''
print(len(nodes))
print(len(shells))
print(len(nodes))
print(len(shells))
'''
Cp = []
with open("Cp.txt","r") as file:
    for line in file:
        Cp.append(float(line))
print(len(Cp))
print(len(shells))

Cp = Cp[:len(mesh.panels)]
print(len(Cp))


V_fs = Vector((1,0,0))
V_fs_norm = V_fs.norm()
mu = []
with open("mu.txt","r") as file:
    for line in file:
        mu.append(float(line))
print(len(mu))
print(len(shells))
mu = mu[:len(shells)]
for j,m in enumerate(mu):
    mesh.panels[j].mu = m
'''
panel = mesh.panels[17]
panel_neighbours = mesh.give_neighbours(panel)

n = len(panel_neighbours)
A = np.zeros((n, 3))
b = np.zeros((n, 1))

for j, neighbour in enumerate(panel_neighbours):
    r_ij = neighbour.r_cp - panel.r_cp
    r_ij = r_ij.transformation(panel.R)
    A[j][0] = r_ij.x
    A[j][1] = r_ij.y
    A[j][2] = r_ij.z
        
    b[j][0] = neighbour.mu - panel.mu
U,S,V = np.linalg.svd(A)
del_mu,_,_,_ = np.linalg.lstsq(A, b, rcond=1e-5)

print(A)
print(U*S*V)
print(b)
print(del_mu)

'''
for i,panel in enumerate(mesh.panels):
            
    panel_neighbours = mesh.give_neighbours(panel)
    panel.sigma = - (panel.n * V_fs)
    panel.Velocity = panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2)
            
    panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2
    if panel.Cp < -0.5:
        print(i)

Cp = [panel.Cp for panel in mesh.panels]
# Surface Contour plot
plot_savedCp_SurfaceContours(mesh.panels, Cp, elevation=20, azimuth=45)
print(len(mesh.panels))