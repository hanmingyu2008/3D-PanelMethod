import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_method_class import PanelMethod, Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from PanelMethod3D.plot_functions import plot_Cp_SurfaceContours
from PanelMethod3D.mesh_class import PanelMesh, Mesh


nodes, shells = 0, 0
nodes = [(-1,-1,-1), (-1,-1,1), (-1,1,-1), (-1,1,1),
         (1,-1,-1), (1,-1,1), (1,1,-1), (1,1,1)]
shells = [[0,1,3,2],[0,4,5,1],[0,2,6,4],[4,6,7,5],[1,5,7,3],[2,3,7,6]]

mesh = PanelMesh(nodes, shells)
mesh = mesh.refinement().refinement().refinement()

V_fs = Vector((1, 0, 0))
panel_method = Steady_Wakeless_PanelMethod(V_fs)
panel_method.solve2(mesh)

q = Vector((0,0,2))
v = Vector((0,0,0))

# print([panel.sigma for panel in mesh.panels])
# print([panel.mu for panel in mesh.panels])
# print([panel.Velocity * panel.n for panel in mesh.panels])
# print([panel.Cp for panel in mesh.panels])

# Surface Contour plot
# plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)