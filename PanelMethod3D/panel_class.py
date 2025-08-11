import numpy as np
from PanelMethod3D.vector_class import Vector
'''
Panel(面元)类和以此为基类的两个类:三角形网格类(triPanel)和四边形网格类(quadPanel)
注意：顶点序列应该是从外部看为逆时针
'''
        
class Panel:
    def __init__(self, position_vector_list):
        self.id = None  # panel的序号
        self.sigma = 0  # panel上的源强度σ,在panel上处处常值
        self.mu = 0  # panel上的偶极子强度σ,在panel上处处常值
        self.num_vertices = len(position_vector_list) # panel顶点个数(只有3或4一般)
        
        # r_vertex[i] : 把panel的各个顶点的坐标排成矩阵一个
        self.r_vertex = np.array(position_vector_list)
        
        self.r_cp = None # panel中心点(控制点)的坐标 
        self.n = None  # panel(中心点处)的单位法向向量
        self.l = None # l和m是(中心点处)的单位切向向量(正交的两个方向)
        self.m = None  # 
        
        # R : 旋转矩阵
        # 下面的两个是把其他点和panel一起旋转了，并且保证panel的中心点r_cp会变成原点
        # r_p_local = R*(r_cpp) = = R*(r_p - r_cp)
        # r_local = R*r_global or r' = R*r  (x' = Ax)
        self.R = None # panel的旋转矩阵，应该是为了方便把panel变成在xy平面上方便计算的
        
        self.r_vertex_local = np.empty_like(self.r_vertex)
        
        self.char_length = 0  # panel的最大对角线距离
        self.area = 0 # panel(表)面积
        
        self.Velocity = Vector((0, 0, 0))  # 中心点处速度矢量
        self.Cp = 0  # 中心点处压力系数Cp Pressure Coefficient
        
        self.set_up_geometry() # 依次调用下面各函数，好在输入panel的顶点列之后可以顺次求出上面除了σ、μ和速度、压力系数之外的大多东西
              
    def set_centroid(self): # 求中心点
        r_cp = Vector((0, 0, 0))
        for i in range(self.num_vertices):
            r_cp = r_cp + self.r_vertex[i]
        
        self.r_cp = r_cp/self.num_vertices
            
    def set_n(self): # 求单位法向量
        r_cp = self.r_cp
        r_vertex = self.r_vertex
        normal = Vector((0, 0, 0))
        for i in range(self.num_vertices):
            j = (i+1)%self.num_vertices
            r_i = r_vertex[i] - r_cp
            r_j = r_vertex[j] - r_cp
            normal = normal + Vector.cross_product(r_i, r_j)
        
        self.area = normal.norm()/2
        self.n = normal/normal.norm()
    
    def set_l(self): # 求单位切向量之一
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_21 = r_1 - r_2
        self.l = r_21/r_21.norm()
    
    def set_m(self): # 求单位切向量之二，其实就是n和l叉乘
        self.m = Vector.cross_product(self.n, self.l)
    
    def set_unit_vectors(self): # 求n l m之综合
        self.set_n()
        self.set_l()
        self.set_m()
            
    def set_R(self): # 旋转矩阵，把l m n分别变成e1 e2 e3或说x' y' z'
        l, m, n = self.l, self.m, self.n
        self.R = np.array([[l.x, l.y, l.z],
                           [m.x, m.y, m.z],
                           [n.x, n.y, n.z]])

    def set_r_vertex_local(self): # 做正交变换之后的panel各个顶点坐标
        n = self.num_vertices
        r_cp = self.r_cp
        r_vertex = self.r_vertex
        r_vertex_local = self.r_vertex_local
        R = self.R
        
        for i in range(n):
            r_vertex_local[i] = r_vertex[i] - r_cp
            r_vertex_local[i] = r_vertex_local[i].transformation(R)

    def set_char_length(self): # 最长对角线(顶点连线)长度
        r = self.r_vertex
        side_lengths = [
            (r[i] - r[i+1]).norm() for i in range(self.num_vertices - 1)
        ]
        
        self.char_length = np.max(side_lengths)      
    
    def set_area(self): # 面积计算 对于平面上多边形A1A2...An,中心点C,相当于求出有向面积S_{CA_iA_{i+1}}求和
        
        r_cp = self.r_cp
        r_vertex = self.r_vertex
        normal = Vector((0, 0, 0))
        for i in range(self.num_vertices):
            j = (i+1)%self.num_vertices
            r_i = r_vertex[i] - r_cp
            r_j = r_vertex[j] - r_vertex
            normal = normal + Vector.cross_product(r_i, r_j)
        
        self.area = normal.norm()/2           
    
    def set_up_geometry(self):
        self.set_centroid()
        self.set_unit_vectors()
        self.set_R()
        self.set_r_vertex_local()
        self.set_char_length()
        # self.set_area()
        pass
   
    
class quadPanel(Panel): # (共面)四边形panel面元
    def __init__(self, vertex0:Vector, vertex1:Vector,
                 vertex2:Vector, vertex3:Vector):
        super().__init__([vertex0, vertex1, vertex2, vertex3])
        
    def set_n(self): # 法向量计算
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_3 = self.r_vertex[2]
        r_31 = r_1 - r_3
        r_4 = self.r_vertex[3]
        r_24 = r_4 - r_2
        cross_product = Vector.cross_product(r_24, r_31)
        n = cross_product/cross_product.norm()        
        self.n = n
        self.area = cross_product.norm()/2
    
    def set_VSAERO_unit_vectors(self):
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_3 = self.r_vertex[2]
        r_4 = self.r_vertex[3]
        r_c = self.r_cp
        
        
        D1 = r_3 - r_1
        D2 = r_4 - r_2
        n = Vector.cross_product(D1, D2)
        self.n = n/n.norm()
        
        m = (r_3 + r_4)/2 - r_c
        self.m = m/m.norm()
        
        self.l = Vector.cross_product(self.m, self.n)
        
        SMP = (r_2 + r_3)/2 - r_c
        self.SMP = SMP.norm() 
        SMQ = (r_3 + r_4)/2 - r_c
        self.SMQ = SMQ.norm()
        
        self.T = (r_3 + r_2)/2 - r_c
    
    def set_char_length(self):
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_3 = self.r_vertex[2]
        r_31 = r_1 - r_3
        r_4 = self.r_vertex[3]
        r_24 = r_4 - r_2
        self.char_length = np.max([r_24.norm(), r_31.norm()])
    
    def set_area(self):    
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_3 = self.r_vertex[2]
        r_4 = self.r_vertex[3]
        r_31 = r_1 - r_3
        r_24 = r_4 - r_2
        cross_product = Vector.cross_product(r_24, r_31)
        self.area = cross_product.norm()/2    
        
    def set_up_geometry(self):
        self.set_centroid()
        self.set_VSAERO_unit_vectors()
        self.set_R()
        self.set_r_vertex_local()
        self.set_char_length()
        self.set_area()    

      
class triPanel(Panel): # 三角形panel面元
    def __init__(self, vertex0:Vector, vertex1:Vector,
                 vertex2:Vector):
        super().__init__([vertex0, vertex1, vertex2])
        
    
    def set_n(self):
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_3 = self.r_vertex[2]
        r_31 = r_1 - r_3
        r_21 = r_1 - r_2
        cross_product = Vector.cross_product(r_21, r_31)
        n = cross_product/cross_product.norm() if cross_product.norm() != 0 else Vector((0,0,0))
        self.n = n
        self.area = cross_product.norm()/2
    
    def set_char_length(self):
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_3 = self.r_vertex[2]
        r_31 = r_1 - r_3
        r_21 = r_1 - r_2
        r_32 = r_2 - r_3
        self.char_length = np.max([r_21.norm(), r_32.norm(), r_31.norm()])
    
    def set_area(self):
            
        r_1 = self.r_vertex[0]
        r_2 = self.r_vertex[1]
        r_3 = self.r_vertex[2]
        r_31 = r_1 - r_3
        r_21 = r_1 - r_2
        cross_product = Vector.cross_product(r_21, r_31)
        self.area = cross_product.norm()/2
    
    def set_up_geometry(self):
        self.set_centroid()
        self.set_unit_vectors()
        self.set_R()
        self.set_r_vertex_local()
        self.set_char_length()
        # self.set_area()
        pass
        
    
if __name__=='__main__':
    from matplotlib import pyplot as plt
    vertex1 = Vector((-1, -1, 2))
    vertex2 = Vector((1, -1, 2))
    vertex3 = Vector((1, 1, -2))
    
    panel = triPanel(vertex1,vertex2,vertex3)
    
    ax = plt.axes(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(-5, 5)
    ax.set_ylim(-5, 5)
    ax.set_zlim(-5, 5)
    # ax.view_init(0, 0)
    x = []
    y = []
    z = []
    for i in range(panel.num_vertices+1):
        i = i % panel.num_vertices
        x.append(panel.r_vertex[i].x)
        y.append(panel.r_vertex[i].y)
        z.append(panel.r_vertex[i].z)
       
    ax.plot3D(x, y, z, color='k', label='panel')
    
    ax.quiver(panel.r_cp.x, panel.r_cp.y, panel.r_cp.z,
              panel.n.x, panel.n.y, panel.n.z,
              color='r', label='normal vector n')
    ax.quiver(panel.r_cp.x, panel.r_cp.y, panel.r_cp.z,
              panel.l.x, panel.l.y, panel.l.z,
              color='b', label='longidutinal vector l')
    ax.quiver(panel.r_cp.x, panel.r_cp.y, panel.r_cp.z,
              panel.m.x, panel.m.y, panel.m.z,
              color='g', label='transverse vector m')
    
    ax.legend()
    
    plt.show()