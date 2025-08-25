import numpy as np
from matplotlib import pyplot as plt

def Wing(filename, length, num, typa = "quad", filled = False):
    # file第一行：名字，譬如 NACA 0012
    # file第二行开始每行两个浮点数作为一个节点的坐标，逆时针排列。另外，第一个和最后一个需要一样！
    coor = []
    lines = []
    with open("coord_seligFmt/"+filename,"r") as filee:
        for line in filee:
            lines.append(line)
    print(lines[0])
    for line in lines[1:]:
        coor.append(tuple([float(x) for x in line.split()]))

    plt.plot([x for x,z in coor],[z for x,z in coor], 'ks-', markersize=1, linewidth=1.5)
    plt.show()
    m = len(coor) - 1
    
    if m%2 != 0:
        raise Exception("我们的2D机翼顶点个数最好是偶数啊(不计重)")
    
    ylist = list(np.linspace(0,length,num+1))
    nodes = []
    for y in ylist:
        for x,z in coor[:-1]:
            nodes.append((x,y,z))
    
    l = m // 2
    x0,z0 = coor[0]
    xk,zk = coor[l]
    semi_xlist = []
    semi_zlist = []
    for i in range(0,l+1):
        x,z = calculate_projection((x0,z0),(xk,zk),coor[i])
        semi_xlist.append(x)
        semi_zlist.append(z)

    for i in range(1,l):
        nodes.append((semi_xlist[i],0,semi_zlist[i]))
    for i in range(1,l):
        nodes.append((semi_xlist[i],length,semi_zlist[i]))
    
    shells = []
    if typa == "quad":
        for i in range(num):
            for k in range(m):
                k_next = (k+1) % m
                shells.append([(i+1)*m+k,(i+1)*m+k_next,i*m+k_next,i*m+k])
    elif typa == "tri":
        for i in range(num):
            for k in range(m):
                k_next = (k+1) % m
                shells.append([(i+1)*m+k,(i+1)*m+k_next,i*m+k_next])
                shells.append([i*m+k,(i+1)*m+k,i*m+k_next])
    else:
        raise Exception("你这个类型我没听过欸")
    
    tot = (num+1)*m

    shells.append([0,1,tot])
    shells.append([0,tot,m-1])
    for i in range(1,l-1):
        shells.append([i,i+1,tot+i,tot+i-1])
        shells.append([m-i-1,m-i,tot+i-1,tot+i])
    shells.append([l-1,l,tot+l-2])
    shells.append([l,l+1,tot+l-2])

    t = num * m
    tot = tot + (l-1)
    shells.append([t+1,t,tot])
    shells.append([t+m-1,tot,t])
    for i in range(1,l-1):
        shells.append([t+i+1,t+i,tot+i-1,tot+i])
        shells.append([t+m-i,t+m-i-1,tot+i,tot+i-1])
    shells.append([t+l,t+l-1,tot+l-2])
    shells.append([t+l+1,t+l,tot+l-2])
    
    mid = num // 2
    shells_to_show = list(range(mid*m,(mid+1)*m))
    if filled:
        return nodes,shells,shells_to_show
    else:
        return nodes[:m*(num+1)],shells[:m*num],shells_to_show
    
def Wing2(filename, length, num, typa = "quad"):
    # 这个和上面那个本质的差别就是不封口，自然而然也不会出现第一个和最后一个是重复的
    coor = []
    lines = []
    with open("coord_seligFmt/"+filename,"r") as filee:
        for line in filee:
            lines.append(line)
    print(lines[0])
    for line in lines[1:]:
        coor.append(tuple([float(x) for x in line.split()]))

    plt.plot([x for x,z in coor],[z for x,z in coor], 'ks-', markersize=1, linewidth=1.5)
    plt.show()
    m = len(coor)
    
    
    ylist = list(np.linspace(0,length,num+1))
    nodes = []
    for y in ylist:
        for x,z in coor:
            nodes.append((x,y,z))
    
    
    shells = []
    if typa == "quad":
        for i in range(num):
            for k in range(m-1):
                k_next = (k+1) 
                shells.append([(i+1)*m+k,(i+1)*m+k_next,i*m+k_next,i*m+k])
    elif typa == "tri":
        for i in range(num):
            for k in range(m-1):
                k_next = (k+1) 
                shells.append([(i+1)*m+k,(i+1)*m+k_next,i*m+k_next])
                shells.append([i*m+k,(i+1)*m+k,i*m+k_next])
    else:
        raise Exception("你这个类型我没听过欸")
    
    
    mid = num // 2
    shells_to_show = list(range(mid*m,(mid+1)*m))
    
    return nodes[:m*(num+1)],shells[:m*num],shells_to_show

def calculate_projection(A, B, C):
    """
    计算点C在由点A和点B形成的直线上的投影坐标
    
    参数:
    A, B, C: 元组或列表，表示平面上的点，格式为 (x, y)
    
    返回:
    投影点的坐标 (x, y)
    """
    # 提取坐标
    Ax, Ay = A
    Bx, By = B
    Cx, Cy = C
    
    # 计算向量AB和AC
    ABx = Bx - Ax
    ABy = By - Ay
    ACx = Cx - Ax
    ACy = Cy - Ay
    
    # 计算点积 AB·AC
    dot_product = ABx * ACx + ABy * ACy
    
    # 计算AB长度的平方
    AB_length_squared = ABx**2 + ABy**2
    
    # 计算投影比例参数t
    if AB_length_squared == 0:  # A和B重合，无法形成直线
        return None
    t = dot_product / AB_length_squared
    
    # 计算投影点坐标
    Px = Ax + t * ABx
    Py = Ay + t * ABy
    
    return (Px, Py)