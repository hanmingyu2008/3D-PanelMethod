import numpy as np
from math import cos,sin,pi
from vector_class import Vector
from panel_method_class import Steady_PanelMethod,VSAERO_panel_velocity,panel_velocity2,panel_velocity3,panel_velocity
from matplotlib import pyplot as plt
from plot_functions import plot_Cp_SurfaceContours,plot_savedCp_SurfaceContours
from mesh_class import PanelMesh
from LSQ import lsq

nodes = [(51.1136,-66.8987,-17.4626),(-52.7887,-64.9883,-17.099),(155.88,-24.6179,-0.00362188),
 (51.9599,-23.5906,-0.00120729),(-51.9599,-22.5633,0.00120729),(-155.344,-19.6652,-10.1751),
 (155.88,23.0396,0.00362188),(51.96,23.0645,0.00120729),(-51.96,23.0894,-0.00120729),(-155.344,23.0676,-10.1732),
 (51.9601,69.7197,0.00362188),(-51.9601,68.7421,-0.00362188)]

nodes2 = [(4103.798, -45.79867245, -235.9695),
(4207.596, -44.813700000000004, -231.024),
(4000.0, 1.02e-05, -250.0),
(4103.798, 0.0491551, -244.868),
(4207.596, 0.0983, -239.736),
(4310.3665, 0.049149575, -224.2925),
(4000.0, 46.783635100000005, -240.91500000000002),
(4103.798, 45.847827550000005, -235.9695),
(4207.596, 44.91202, -231.024),
(4310.366499999999, 41.9974022875, -216.14175),
(4103.798, 91.6465, -227.07100000000003),
(4207.596, 89.72574, -222.312)]
shells = [[4,1,0,3],[2,6,7,3],[4,3,7,8],[4,8,9,5],[11,8,7,10]]

mesh = PanelMesh(nodes,shells)

mu = [25.425937892079475, 32.87213499130843, 25.42983296682372, 12.62608022674372, 25.488996085081688]

panel = mesh.panels[2]
print(panel.n)

V_fs = Vector((cos(pi/60),0,sin(pi/60)))
V_fs_norm = V_fs.norm()

for i in range(len(mesh.panels)):
    mesh.panels[i].mu = mu[i]

panel_neighbours = mesh.give_neighbours3(panel, 3)
panel.Velocity = panel_velocity3(panel, panel_neighbours, V_fs, rcond=1e-2)
print(1-(panel.Velocity.norm())**2)

panel_neighbours = mesh.give_neighbours(panel)
v = VSAERO_panel_velocity(V_fs, panel, panel_neighbours)
print(v)
print(1-(v.norm())**2)

v2 = panel_velocity(panel,mesh.give_neighbours(panel),V_fs)
print(v2)
print(1-(v2.norm())**2)
plot_savedCp_SurfaceContours(mesh.panels, mu)