import numpy as np
from vector_class import Vector
from influence_coefficient_functions import influence_coeff
from influence_coefficient_functions import compute_source_panel_velocity,compute_dipole_panel_velocity
from mesh_class import PanelMesh
from panel_class import quadPanel, triPanel, Panel

## Lothar的论文里面有几个小测试，通过比对数据我觉得我前面的代码可以通过考量

panel1 = quadPanel(Vector((-0.5,-0.5,0)),Vector((0.5,-0.5,0)),Vector((0.5,0.5,0)),Vector((-0.5,0.5,0)))
panel2 = quadPanel(Vector((-0.7,-1,0)),Vector((-1,1,0)),Vector((1,0.5,0.8)),Vector((1,-1,1)))
panel3 = triPanel(Vector((0,0,0)),Vector((1,0,0)),Vector((1,1,0)))

q11 = Vector((0,0,2))
q12 = Vector((0,0,-0.5))
q13 = Vector((0.5,1,0))
q14 = Vector((-0.25,0.25,0.5))
q15 = Vector((0.1,0.4,8))

q1 = [None,q11,q12,q13,q14,q15]
'''
vs, vd = compute_source_panel_velocity(q1[3],panel1,1), compute_dipole_panel_velocity(q1[3],panel1,1)
ps, pd = influence_coeff(q1[3],panel1)

print(ps,vs)
print(pd,vd)
'''
q21 = Vector((0,0,2))
q22 = Vector((1,-1,0.812))
q23 = Vector((2,-0.11,0.1))
q24 = Vector((0.5,-0.5,-0.1))
q25 = Vector((1.5,-0.5,-9))

q2 = [None,q21,q22,q23,q24,q25]

i=5

vs, vd = compute_source_panel_velocity(q2[i],panel2,1), compute_dipole_panel_velocity(q2[i],panel2,1)
ps, pd = influence_coeff(q2[i],panel2)

q31 = Vector((0,0,2))
q32 = Vector((0,0,-0.5))
q33 = Vector((0.5,1,0))
q34 = Vector((-0.25,0.25,0.5))
q35 = Vector((0.1,0.4,8))

q3 = [None,q31,q32,q33,q34,q35]

for i in range(1,6):

    vs, vd = compute_source_panel_velocity(q3[i],panel3,1), compute_dipole_panel_velocity(q3[i],panel3,1)
    ps, pd = influence_coeff(q3[i],panel3)

    print(ps,vs)
    print(pd,vd)