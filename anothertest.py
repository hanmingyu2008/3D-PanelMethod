import numpy as np
from math import cos,sin,pi
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_method_class import Steady_PanelMethod,VSAERO_panel_velocity,panel_velocity3
from matplotlib import pyplot as plt
from PanelMethod3D.plot_functions import plot_Cp_SurfaceContours,plot_savedCp_SurfaceContours
from PanelMethod3D.mesh_class import PanelAeroMesh
from PanelMethod3D.LSQ import lsq


nodes = []
shells = []
shells_ids = {"body":[],"wake":[]}

with open("plane_refined.txt","r") as filee:
    for line in filee:
        if line[0] == " ":
            x,y,z = [float(t) for t in line.split()]
            nodes.append((x,y,z))
        if line[0] == "1" and line[1] != "0":
            labels = [int(t) for t in line[1:].split()]
            shells_ids["body"].append(len(shells))
            shells.append(labels[:4])
        if line[0] == "1" and line[1] == "0":
            labels = [int(t) for t in line[2:].split()]
            shells_ids["wake"].append(len(shells))
            shells.append(labels[:4])
        if line[0] == "2":
            labels = [int(t) for t in line[1:].split()]
            shells_ids["body"].append(len(shells))
            shells.append(labels[:3])

mesh = PanelAeroMesh(nodes, shells, shells_ids)

V_fs = Vector((cos(pi/60),0,sin(pi/60)))
V_fs_norm = V_fs.norm()
'''
mesh = mesh.refinement()

with open("plane_refined.txt","w") as filee:
    for x,y,z in mesh.nodes:
        print(" ",x,y,z,file=filee)
    for i in mesh.shells_ids["body"]:
        lst = mesh.shells[i]
        if len(lst) == 3:
            a,b,c = lst
            print(2,a,b,c,file=filee)
        if len(lst) == 4:
            a,b,c,d = lst
            print(1,a,b,c,d,file=filee)
    for i in mesh.shells_ids["wake"]:
        lst = mesh.shells[i]
        a,b,c,d = lst
        print(10,a,b,c,d,file=filee)
'''       

'''
method = Steady_PanelMethod(V_fs)
method.solve(mesh)

with open("mu.txt","w") as filee:
    for panel in mesh.panels:
        print(panel.mu,file=filee)
with open("Cp.txt","w") as filee:
    for panel in mesh.panels:
        print(panel.Cp,file=filee)
'''
Cp = []
mu = []
with open("Cp_new.txt","r") as filee:
    for line in filee:
        Cp.append(float(line))
with open("mu_refine.txt","r") as filee:
    for line in filee:
        mu.append(float(line))
L = 0

for i,panel in enumerate(mesh.panels):
    panel.mu = mu[i]
    panel.Cp = Cp[i]

# plot_savedCp_SurfaceContours(mesh.panels[:len(mesh.shells_ids["body"])],Cp[:len(mesh.shells_ids["body"])])


lst = mesh.give_neighbours(mesh.panels[8227])
lst.append(mesh.panels[8227])
plot_savedCp_SurfaceContours(lst,[panel.mu for panel in lst])
    

for i in mesh.shells_ids["body"]:
    panel = mesh.panels[i]
    L -= Cp[i] * panel.area * (panel.n * Vector((0,0,1)))
print(L)


L = 0
L_ref = 0
for i in [8227]: # mesh.shells_ids["body"]:
    panel = mesh.panels[i]
    panel.sigma = -(panel.n * V_fs)
    neighbours = mesh.give_neighbours(panel)
    if len(neighbours) == 4:
        b = True # 用来判断到底可不可以用VSAERO来计算
        ok = [True]*4
        for j,neigh in enumerate(neighbours):
            if neigh.n * panel.n < 0.9:
                ok[j] = False
                if neighbours[(j+2)%4].num_vertices == 3:
                    b = False
                else:
                    nei_new = mesh.neighbour_of_neighbour(panel,(j+2)%4)
                    if nei_new.n * panel.n < 0.9:
                        b = False
                    else:
                        neighbours[j] = nei_new
        if not ((ok[0] or ok[2]) and (ok[1] or ok[3])):
            b = False

        if b:
            print(i)
            v = panel_velocity3(panel, mesh.give_neighbours3(panel,3), V_fs, rcond = 1e-2)
            panel.Velocity = VSAERO_panel_velocity(V_fs,panel,neighbours,ok[0],ok[1],ok[2],ok[3])
            print(v.transformation(panel.R))
            print(panel.Velocity.transformation(panel.R))
            cp = 1 - (v.norm()/V_fs_norm)**2
            Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2
            print(abs(cp-Cp),abs(cp-Cp)/abs(cp))
        else:
            panel_neighbours = mesh.give_neighbours3(panel,3)
            panel.Velocity = panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2)
    else:
            panel_neighbours = mesh.give_neighbours3(panel,3)
            panel.Velocity = panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2)        
            
    panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2

for i in mesh.shells_ids["body"]:
    panel = mesh.panels[i]
    L -= panel.Cp * panel.area * (panel.n * Vector((0,0,1)))
    L_ref += panel.area * abs(panel.n * Vector((0,0,1)))

print(L)
print(L_ref)