'''
对于一个.ply的表面网格文件,当我们求出Cp值之后可以把它作为.ply文件里面的face quality写进去,以便在MeshLab中可视化观察
'''

filename = "plane_0001_cfd"
lines = []
with open(filename+".ply","r") as filee:
    for line in filee:
        lines.append(line)
ind = lines.index("end_header\n")
a = 10386 
b = 20632

Cp = []
with open(filename+"_Cp.txt","r") as filee:
    for line in filee:
        Cp.append(float(line))

with open(filename+"_withCp.ply","w") as filee:
    for i in range(ind):
        print(lines[i][:-1],file=filee)
    print("property float quality",file=filee)
    print(lines[ind][:-1],file=filee)
    for i in range(a):
        print(lines[ind+1+i][:-1],file=filee)
    for i in range(b):
        print(lines[ind+a+1+i][:-1],Cp[i],file=filee)