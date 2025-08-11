import numpy as np
import matplotlib.pyplot as plt
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.influence_coefficient_functions import influence_coeff
from PanelMethod3D.influence_coefficient_functions import compute_source_panel_velocity,compute_dipole_panel_velocity
from PanelMethod3D.mesh_class import PanelMesh,PanelAeroMesh
from PanelMethod3D.panel_method_class import Steady_Wakeless_PanelMethod,Steady_PanelMethod
from PanelMethod3D.plot_functions import plot_Cp_SurfaceContours
from PanelMethod3D.influence_coefficient_functions import influence_coeff

epsilon = 0
nodes = [(0,0,1),(-1,0,0),(0,-1,0),(1,0,0),(0,1,0),(0,0,-1),(-1,-1,epsilon),(-1,-1,-epsilon)]
shells = [[0,1,2],[0,2,3],[0,3,4],[0,4,1],[2,1,5],[3,2,5],[4,3,5],[1,4,5],[1,2,6],[1,6,2]]#,[1,6,7],[2,7,6]]

mesh = PanelMesh(nodes,shells)
method = Steady_Wakeless_PanelMethod(Vector((1,0,4)))

# method.solve(mesh)

print(influence_coeff(mesh.panels[0].r_cp,mesh.panels[-2]))
print(influence_coeff(mesh.panels[0].r_cp,mesh.panels[-1]))

for panel in mesh.panels:
    panel.Cp = panel.mu
plot_Cp_SurfaceContours(mesh.panels)