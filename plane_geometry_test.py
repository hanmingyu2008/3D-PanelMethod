import numpy as np
from vector_class import Vector
from panel_method_class import panel_velocity2,panel_velocity,panel_velocity3,panel_velocity_new,VSAERO_panel_velocity
from matplotlib import pyplot as plt
from plot_functions import plot_savedCp_SurfaceContours,plot_Cp_SurfaceContours
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
        if line[0] == "1" and line[1] == "0":
            labels = [int(t)-1 if t not in ["2071","2158"] else 1259 for t in line[2:].split()]
            shells_ids["wake"].append(len(shells))
            shells.append(labels[:4])
        if line[0] == "2":
            labels = [int(t)-1 if t not in ["2071","2158"] else 1259 for t in line[1:].split()]
            shells_ids["body"].append(len(shells))
            shells.append(labels[:3])

mesh = PanelAeroMesh(nodes, shells, shells_ids)

mesh2 = mesh.refinement()
for panel_id in mesh2.TrailingEdge["suction side"]:
    mesh2.panels[panel_id].Cp = 1
print(mesh2.wake_sheddingPanels)
print(len(mesh2.shells_ids["body"]))
print(len(shells_ids["wake"]))

# plot_Cp_SurfaceContours(mesh2.panels)