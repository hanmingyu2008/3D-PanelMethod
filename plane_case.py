import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_method_class import Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from PanelMethod3D.mesh_class import PanelMesh
import time

## READ FROM OFF FILES
filename = "airplane_0032_nowheel"
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