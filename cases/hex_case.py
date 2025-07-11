import numpy as np
from vector_class import Vector
from panel_method_class import Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from plot_functions import plot_Cp_SurfaceContours
from mesh_class import PanelMesh

radius = 1
num_longitude, num_latitude = 21, 20
nodes = [(1,1,1),
         (1,1,-1),
         (1,-1,1),
         (1,-1,-1),
         (-1,1,1),
         (-1,1,-1),
         (-1,-1,1),
         (-1,-1,-1)]
shells = [[0,1,3],[0,2,3],
          [0,2,6],[0,4,6],
          [0,1,5],[0,4,5],
          [1,3,5],[3,5,7],
          [2,3,6],[3,6,7],
          [4,5,6],[5,6,7]]

mesh = PanelMesh(nodes, shells)

# mesh = PanelMesh.generate_from_stl_file("Sphere_0080")


V_fs = Vector((1, 0, 0))

panel_method = Steady_Wakeless_PanelMethod(V_fs)
panel_method.solve(mesh)

# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)