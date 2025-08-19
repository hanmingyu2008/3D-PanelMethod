import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.influence_coefficient_functions import influence_coeff
from PanelMethod3D.influence_coefficient_functions import compute_source_panel_velocity,compute_dipole_panel_velocity
from PanelMethod3D.mesh_class import PanelMesh
from PanelMethod3D.panel_class import quadPanel, triPanel, Panel

panel = quadPanel(Vector((-1,-1,0)),Vector((-1,1,0)),Vector((1,1,0)),Vector((1,-1,0)))
panel1 = quadPanel(Vector((-1,-1,0)),Vector((-1,0,0)),Vector((0,0,0)),Vector((0,-1,0)))
panel2 = quadPanel(Vector((-1,1,0)),Vector((0,1,0)),Vector((0,0,0)),Vector((-1,0,0)))
panel3 = quadPanel(Vector((1,1,0)),Vector((1,0,0)),Vector((0,0,0)),Vector((0,1,0)))
panel4 = quadPanel(Vector((1,-1,0)),Vector((0,-1,0)),Vector((0,0,0)),Vector((1,0,0)))

q = Vector((0.5,0.5,1))
vs = compute_source_panel_velocity(q,panel,1.0)
vssum = compute_source_panel_velocity(q,panel1,1.0)+compute_source_panel_velocity(q,panel2,1.0)+\
    compute_source_panel_velocity(q,panel3,1.0)+compute_source_panel_velocity(q,panel4,1.0)

print(vs)
print(vssum)

vd = compute_dipole_panel_velocity(q,panel,1.0)
vdsum = compute_dipole_panel_velocity(q,panel1,1.0)+compute_dipole_panel_velocity(q,panel2,1.0)+\
    compute_dipole_panel_velocity(q,panel3,1.0)+compute_dipole_panel_velocity(q,panel4,1.0)

print(vd)
print(vdsum)