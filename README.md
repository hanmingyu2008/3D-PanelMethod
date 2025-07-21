# 3D-PanelMethod
A panel method solver for 3d incompressible inviscid potential flow.
Wake haven't been included in consideration, but I will try to add it later.

## Explaination of Algorithm Employed

The method we employ can be find in Chapter 3,9,12 of 《Low Speed Aerodynamics: From Wing Theory To Panel Method》 written by Joseph Katz. In calculation of the matrixs B,C, we used formula from 《A Comprehensive and Practical Guide to the Hess and Smith Constant Source and Dipole Panel》 written by Lothar Birk.

The Algorithm can be described precisely as follow:

First we assume there is a source of strength $\sigma$ and a doublet of strength $\mu$ on every point of the surface, and $\sigma,\mu$ are constant on each panel. The equation we have is that normal velocity on every point on the surface is 0. That is:
$$\frac{\partial}{\partial\eta}\Phi = \vec{v}\cdot\eta = 0 \forall x\in S$$

After process in Katz., we can reach an algebraic system of $\vec{\sigma},\vec{\mu}$ (vector composed of $\sigma$ and $\mu$ in every panel):
$$B\vec{\sigma}+C\vec{\mu}=0$$
$$\sigma_i = \vec{V}_\infty \cdot \eta_i$$
Here $\sigma_i$ and $\eta_i$ denotes $\sigma$ and outer normal vector on panel-i.

So $\sigma,\mu$ on every panel can be calculated, and then we can see (again, see more detailed explaination in Katz's book)
$$\vec{l}\cdot\triangledown\Phi = -\frac{\partial\mu}{\partial\vec{l}}$$
for every direction $\vec{l}$ especially those tangential to the surface. Then the problem of calculating Cp is turned into directional deriatives on surface, which is always solved by LeastSquare solver (see LSQ.py). (Remeber that the velocity on normal direction is always zero!)

Once we get velocity control point of every panels, we can directly calculate pressure coefficient Cp. And therefore, we can get force on every direction($\vec{i}$) through:
$$F_i = -\sum\limits_{k:panel} Cp_kS_k(\eta_k\cdot\vec{i})q$$ 

So the whole algorithm can be splitted into several parts:
step 0. Construct the mesh, written it into forms of vortexs and panels.
step 1. Compose the matrixs B,C in the algebraic system using formula provided by Mr. Lothar Birk.
step 2. Solve the equation to get $\mu_i,\sigma_i$
step 3. using $\sigma,\mu$ we just get to calculate velocity on control point of each panels,and therefore reach Cp on each panels.
step 4. Calculate forces such as lifting forces. 

## Structure
Now we explain the structure of our code in this repository.

FOR STEP ZERO 

vector_class.py: A class Vector and simple implementation of algebraic calculation.
panel_class.py: A class Panel and two class "triPanel","quadPanel" based on class "Panel".
mesh_class.py: A class Mesh, composed of several(maybe thousands of) Panels

FOR STEP ONE

influence_coefficient_functions.py: Implemented formulas in Lothar's paper, being an auxillary function for the construction of B,C. We use "is_inside_polygon.py" here.
panel_method_class.py: A class of PanelMethod, static method "influence_coeff_matrices" composes B,C.

FOR STEP TWO

panel_method_class.py: In function solve, or functions with similar name such as solve2.

FOR STEP THREE

panel_method_class.py: In function panel_velocity, or functions with similar name such as panel_velocity2.

FOR STEP FOUR

The calculation of lifting force isn't implemented.

OTHER PARTS

mu_plot.py: to show the distribution of mu on surface, after Step Two.
apame.py: used in function panel_velocity_apame of PanelMethod Class.
plot_functions.py: functions to plot the distribution of Cp using matplotlib.
LSQ.py: solving LeastSquare Problem.

SOME TESTS:

1. SPHERE
We designed a sphere class to more precisely construct panels and mesh on a sphere.
And solve the invisid potential problem on the surface of a sphere in sphere_case.py. Also compare it with the analytic solution.

2. HEX
A simple hex to solve.

3. PLANE
This is the plane example on 3dpanelmethod.com.(The mesh isn't defined so well, but the result can be compared with pictures on this website).
wake_test_fuselage_wake.inp is the coordinates of each vortex and vortex_ids of each panel.
In plane case, we assumed you have calculated mu out through other computers, so we provided a mu.txt, and plane_case.py simply implemented steps after the second.

More details could be found in python codes.

Thanks for help from code of Mr. @ichamial