import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_method_class import Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from PanelMethod3D.mesh_class import PanelMesh
import time
from PanelMethod3D.plot_functions import plot_savedCp_SurfaceContours
## READ FROM OFF FILES
filename = "1"
filetype = "off"
with open(filename+"."+filetype,"r") as filee:
    lines =[]
    for line in filee:
        lines.append(line)
nodes = []
shells = []
a,b,c = [int(x) for x in lines[1].split()]
for ind in range(2,a+2):
    x,y,z = [float(x) for x in lines[ind].split()]
    nodes.append((x,y,z))
for ind in range(a+2,a+b+2):
    shells.append([int(x) for x in lines[ind].split()[1:]])
print("End reading",flush=True)
time1 = time.time()
mesh = PanelMesh(nodes,shells)
time2 = time.time()
print(f"End Construction Meshs, time cost: {time2-time1:.2f}",flush=True)

method = Steady_Wakeless_PanelMethod(Vector((1,0,0)))
method.solve2(mesh)
with open(filename+"_Cp.txt","w") as filee:
    for panel in mesh.panels:
        print(panel.Cp,file=filee)
with open(filename+"_mu.txt","w") as filee:
    for panel in mesh.panels:
        print(panel.mu,file=filee)
'''
mu = []
with open(filename+"_mu.txt","r") as filee:
    for line in filee:
        mu.append(float(line))
Cp = []
with open(filename+"_Cp.txt","r") as filee:
    for line in filee:
        Cp.append(float(line))
plot_savedCp_SurfaceContours(mesh.panels,mu)
'''