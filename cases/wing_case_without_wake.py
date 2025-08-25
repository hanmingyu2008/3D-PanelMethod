from PanelMethod3D.mesh_class import PanelMesh
from wing_class import Wing
from PanelMethod3D.panel_method_class import Steady_Wakeless_PanelMethod
from PanelMethod3D.vector_class import Vector
from matplotlib import pyplot as plt
from PanelMethod3D.plot_functions import plot_Cp_SurfaceContours


nodes,shells,shells_to_show = Wing("naca0012_new.dat",1,10)
print(len(nodes),len(shells))
mesh = PanelMesh(nodes,shells)


V_fs = Vector((1, 0, 0))
panel_method = Steady_Wakeless_PanelMethod(V_fs)
panel_method.solve_newvelo(mesh)

# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=-150, azimuth=-120)

xlist = []
Cplist = []


for id in shells_to_show:
    xlist.append(mesh.panels[id].r_cp.x)
    Cplist.append(mesh.panels[id].Cp)

plt.plot(xlist, Cplist, 'ks--', markerfacecolor='r', label='Panel Method')
plt.legend()
plt.grid()
plt.gca().invert_yaxis()
plt.show()
plt.show()