import numpy as np
from vector_class import Vector
from panel_method_class import PanelMethod, Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from plot_functions import plot_Cp_SurfaceContours
from mesh_class import PanelMesh, Mesh
from sphere import sphere


radius = 1
num_longitude, num_latitude = 21, 20
nodes, shells = sphere(radius, num_longitude, num_latitude,
                                    mesh_shell_type='quadrilateral')

mesh = PanelMesh(nodes, shells)

# mesh = PanelMesh.generate_from_stl_file("Sphere_0080")

V_fs = Vector((1,0,0))
method = Steady_Wakeless_PanelMethod(V_fs)

mesh = PanelMesh(nodes, shells)

method.solve(mesh)

Cp = [panel.Cp for panel in mesh.panels]
# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)