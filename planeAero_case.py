import numpy as np
from vector_class import Vector
from panel_method_class import panel_velocity2,panel_velocity,panel_velocity3,panel_velocity_new,VSAERO_panel_velocity
from matplotlib import pyplot as plt
from plot_functions import plot_savedCp_SurfaceContours
from mesh_class import PanelAeroMesh
from LSQ import lsq


nodes = []
shells = []
shells_ids = {"body":[],"wake":[]}

with open("wake_test_fuselage_wake.inp","r") as filee:
    for line in filee:
        if line[0] == " ":
            x,y,z = [float(t) for t in line.split()]
            nodes.append((x,y,z))
        if line[0] == "1" and line[1] != "0":
            labels = [int(t)-1 if t not in ["2071","2158"] else 1259 for t in line[1:].split()]
            shells_ids["body"].append(len(shells))
            shells.append(labels[:4])
        '''
        if line[0] == "1" and line[1] == "0":
            labels = [int(t)-1 if t not in ["2071","2158"] else 1259 for t in line[2:].split()]
            shells_ids["wake"].append(len(shells))
            shells.append(labels[:4])'''
        if line[0] == "2":
            labels = [int(t)-1 if t not in ["2071","2158"] else 1259 for t in line[1:].split()]
            shells_ids["body"].append(len(shells))
            shells.append(labels[:3])

mesh = PanelAeroMesh(nodes, shells, shells_ids)

V_fs = Vector((1,0,0))
V_fs_norm = V_fs.norm()
mu = []         ## 先在集群上计算出mu的值再回来计算Cp
Cp = []                  
with open("mu.txt","r") as file:
    for line in file:
        mu.append(float(line))
with open("Cp.txt","r") as file:
    for line in file:
        Cp.append(float(line))
print(len(mu))
print(len(shells))
mu = mu[:len(shells)]
for j,m in enumerate(mu):
    mesh.panels[j].mu = m
L = 0
L_ref = 0
for i,panel in enumerate(mesh.panels):
    if i == 1905:
        pass
    panel.sigma = -(panel.n * V_fs)
    neighbours = mesh.give_neighbours(panel)
    if len(neighbours) == 4:
        b = True
        for neigh in neighbours:
            if neigh.num_vertices != 4:
                b = False
            if neigh.n * panel.n < 0.7:
                    b = False
        if b:
            v = VSAERO_panel_velocity(V_fs,panel,neighbours)
            panel_neighbours = mesh.give_neighbours3(panel,3)
            panel.Velocity = v # panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2) 
            cp = 1 - (v.norm()/V_fs_norm)**2
            Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2
            # print(i,abs(cp-Cp),abs(cp-Cp)/abs(Cp))
        else:

            panel_neighbours = mesh.give_neighbours3(panel,3)
            panel.Velocity = panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2)
    else:
            panel_neighbours = mesh.give_neighbours3(panel,3)
            panel.Velocity = panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2)        
            
    panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2

# VSAERO_onbody_analysis(V_fs, mesh)
Cp = [panel.Cp for panel in mesh.panels]

for panel in mesh.panels:
    L -= panel.Cp * panel.area * (panel.n * Vector((0,0,1)))
    L_ref += panel.area * abs(panel.n * Vector((0,0,1)))

# 335的邻居:1832，336，353，334
panels = [mesh.panels[1905],mesh.panels[1865],mesh.panels[1906],mesh.panels[1923],mesh.panels[1904]]
plot_savedCp_SurfaceContours(mesh.panels, Cp, elevation=20, azimuth=45)
print(L)
print(L_ref)
