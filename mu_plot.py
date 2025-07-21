from plot_functions import plot_Cp_SurfaceContours
from mesh_class import PanelMesh

# 其实很简单,就是把panel的Cp属性赋以mu的值，这样原先用来画出Cp图像的函数就可以用来画出mu的图像。
# mu应当是一个比较平滑的函数，如果mu不太平滑，那么可能完蛋了。
nodes = []
shells = []

with open("wake_test_fuselage_wake.inp","r") as filee:
    for line in filee:
        if line[0] == " ":
            x,y,z = [float(t) for t in line.split()]
            nodes.append((x,y,z))
        if line[0] == "1" and line[1] != "0":
            labels = [int(t)-1 for t in line[1:].split()]
            shells.append(labels[:4])
        if line[0] == "2":
            labels = [int(t)-1 for t in line[1:].split()]
            shells.append(labels[:3])
mesh = PanelMesh(nodes, shells)


mu = []
with open("mu.txt","r") as file:
    for line in file:
        mu.append(float(line))
print(len(mu))
print(len(shells))

mu = mu[:len(mesh.panels)]
print(len(mu))

for i,panel in enumerate(mesh.panels):
            
    panel.Cp = mu[i]

# Surface Contour plot
plot_Cp_SurfaceContours(mesh.panels, elevation=20, azimuth=45)
print(len(mesh.panels))