import numpy as np
from vector_class import Vector
from panel_method_class import PanelMethod, Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from plot_functions import plot_savedCp_SurfaceContours
from mesh_class import PanelMesh, Mesh
from sphere import sphere


radius = 1
num_longitude, num_latitude = 61, 60
nodes, shells = sphere(radius, num_longitude, num_latitude,
                                    mesh_shell_type='quadrilateral')

mesh = PanelMesh(nodes, shells)

# mesh = PanelMesh.generate_from_stl_file("Sphere_0080")


print(len(nodes))
print(len(shells))

Cp = []
with open("Cp.txt","r") as file:
    for line in file:
        Cp.append(float(line))
print(len(Cp))

mesh = PanelMesh(nodes, shells)

# Surface Contour plot
plot_savedCp_SurfaceContours(mesh.panels, Cp, elevation=20, azimuth=45)