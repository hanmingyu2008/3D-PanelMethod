import numpy as np
import matplotlib.pyplot as plt
from vector_class import Vector
from influence_coefficient_functions import influence_coeff
from influence_coefficient_functions import compute_source_panel_velocity,compute_dipole_panel_velocity
from mesh_class import PanelMesh,PanelAeroMesh
from panel_method_class import Steady_Wakeless_PanelMethod,Steady_PanelMethod
from plot_functions import plot_Cp_SurfaceContours
from sphere import sphere


nodes = [(0,0,0),(0,0,1),(0,0,-1),(-1,0,0),(-1,0,1),(-1,0,-1),(1,0,0),(1,0,1),(1,0,-1),(0,-1,0),(0,-2,0),
         (-1,-1,0),(-1,-2,0),(1,-1,0),(1,-2,0),(0,-3,0),(-1,-3,0),(1,-3,0)]
shells = [[0,1,4,3],[0,6,7,1],[0,2,8,6],[0,3,5,2],[0,3,11,9],[6,0,9,13],[9,11,12,10],[13,9,10,14],
          [10,12,16,15],[14,10,15,17]]
shells_ids = {"body":[0,1,2,3],"wake":[4,5,6,7,8,9]}
mesh = PanelAeroMesh(nodes,shells,shells_ids)

for id in mesh.TrailingEdge["pressure side"]:
    mesh.panels[id].Cp = 1

print(mesh.wake_sheddingPanels)
print(mesh.shell_neighbours)
# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)