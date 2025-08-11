import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_method_class import Steady_Wakeless_PanelMethod
from matplotlib import pyplot as plt
from PanelMethod3D.plot_functions import plot_Cp_SurfaceContours
from PanelMethod3D.mesh_class import PanelMesh
from cases.sphere import sphere


radius = 1
num_longitude, num_latitude = 20, 21
nodes, shells = sphere(radius, num_longitude, num_latitude,
                                    mesh_shell_type='quadrilateral')

mesh = PanelMesh(nodes, shells)


V_fs = Vector((0,0,1))
panel_method = Steady_Wakeless_PanelMethod(V_fs)
panel_method.solve_newvelo(mesh)


saved_ids = []
for panel in mesh.panels:
    if abs(panel.r_cp.y) <= 10**(-3):
        saved_ids.append(panel.id)

r = 1 
x0, y0 = 0, 0 
analytical_theta = np.linspace(-np.pi, np.pi, 200)
analytical_cp = 1 - (3/2*np.sin(analytical_theta))**2
fig = plt.figure()
plt.plot(analytical_theta*(180/np.pi), analytical_cp ,'b-',
            label='Analytical - sphere')

thetas = []
Cp = []
for id in saved_ids:
    # print(mesh.panels[id].r_cp)
    theta = np.arctan2(mesh.panels[id].r_cp.x, mesh.panels[id].r_cp.z)
    thetas.append(np.rad2deg(theta))
    Cp.append(mesh.panels[id].Cp)
    
plt.plot(thetas, Cp, 'ks', markerfacecolor='r',
            label='Panel Method - Sphere')
plt.xlabel("angle [degrees]")
plt.ylabel("Cp")
plt.title("Cp distribution along the great circle of a sphere")
plt.legend()
plt.grid()
plt.show()

# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)