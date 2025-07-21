import numpy as np
from vector_class import Vector
from influence_coefficient_functions import influence_coeff
from influence_coefficient_functions import compute_source_panel_velocity,compute_dipole_panel_velocity
from mesh_class import PanelMesh
from LSQ import LeastSquares,lsq
from apame import cal_velocity
## 最关键的函数了，这里面实现了多种用于求解势流问题的方法，包括solve和solve2等等
# 这些方法各有千秋，其实并不算很难看懂，但我们也对每个函数都进行解释。
# 部分solve们是参照README里面说的方法进行的，一直到计算出mu都是一摸一样，粘贴复制的。所以我们不再介绍这一部分了。
# 关键在于如何从mu到速度
# 也有少部分sovle采用的是别的方法，这些应该改都标注了"这个函数尚未改对"之类的标注，我们也会再详细解释。

class PanelMethod:
    
    def __init__(self, V_freestream):
        self.V_fs = V_freestream  
    
    def set_V_fs(self, Velocity, AngleOfAttack, SideslipAngle):
        alpha = np.deg2rad(AngleOfAttack)
        beta = np.deg2rad(SideslipAngle)
        Vx = Velocity * np.cos(alpha) * np.cos(beta)
        Vy = Velocity * np.cos(alpha) * np.sin(beta)
        Vz = - Velocity * np.sin(alpha)
        self.V_fs = Vector((Vx, Vy, Vz))


class Steady_Wakeless_PanelMethod(PanelMethod): # 或许什么时候我们会把Wake给补上，现在只好先叫做Wakeless
    
    def __init__(self, V_freestream):
        super().__init__(V_freestream)
    
    def solve(self, mesh:PanelMesh):
        
        for panel in mesh.panels:
            panel.sigma = source_strength(panel, self.V_fs) 
        
        B, C = self.influence_coeff_matrices(mesh.panels)
        
        RHS = right_hand_side(mesh.panels, B)
        
        doublet_strengths = np.linalg.solve(C, RHS)
        
        for panel_id, panel in enumerate(mesh.panels):
            panel.mu = doublet_strengths[panel_id]
        
        V_fs_norm = self.V_fs.norm()
        
        for panel in mesh.panels:
            ## 速度计算的方法
            # 把panel的所有公共边邻居都取出来，三个维度上列出来最小二乘法，采用稍微不太科学的np自带函数求解。
            # 什么叫"三个维度上最小二乘法"见函数panel_velocity2
            # 其实这个方法还不错，可以解决很多情况，只会在边角处出现点问题吧。
            # 其实吧panel_velocity2换成panel_velocity也还可以，也是会在边角处有点小问题。
            # 要我说有点小问题(在小panel上有点小偏差不是大问题嘿)
            panel_neighbours = mesh.give_neighbours(panel)
            panel.Velocity = panel_velocity2(panel, panel_neighbours, self.V_fs, rcond = 1e-3)
            
            panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2 

    def solve_new(self, mesh:PanelMesh):  # 这个函数还没有改对！！！不要使用！！！
        
        for panel in mesh.panels:
            panel.sigma = source_strength(panel, self.V_fs) 
        
        B, C = self.influence_coeff_matrices(mesh.panels)
        
        RHS = right_hand_side(mesh.panels, B)
        
        doublet_strengths = np.linalg.solve(C, RHS)
        
        for panel_id, panel in enumerate(mesh.panels):
            panel.mu = doublet_strengths[panel_id]
        
        V_fs_norm = self.V_fs.norm()
        
        for panel in mesh.panels:
            ## 速度计算方法
            # 把我们所有的panel对这个控制点的速度贡献全都加起来不就好了吗
            # 它就是不知道问题出现在哪里啊
            # 属于是前面根本就没有变化，我就改了这一步，他竟然就不对了??
            # 别说Sphere了，hex的问题，都没有办法使得panel.Velocity * panel.n = 0
            panel.Velocity = panel_velocity_new(panel, mesh, self.V_fs)            
            
            panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2

    def solve_newvelo(self, mesh:PanelMesh): # 这个函数还没有改对！！！不要使用！！！
        # 这个函数的B,C和其他的几个并不一样，完全仿照二维的panel method进行的
        # 但是他们不这么干可能也是有他们的道理，这个东西就不指望改对了吧
        # 要说战绩，是能保证panel.Velocity * panel.n = 0,毕竟我方程列的就是这个
        # 但是sphere和解析解对不上啊
        for panel in mesh.panels:
            panel.sigma = source_strength(panel, self.V_fs) 
        
        B, C = self.influence_velocoeff_matrices(mesh.panels)
        
        RHS = right_hand_side(mesh.panels, B)

        for panel_id, panel in enumerate(mesh.panels):
            RHS[panel_id] -= panel.n * self.V_fs
        
        doublet_strengths = np.linalg.solve(C, RHS)
        
        for panel_id, panel in enumerate(mesh.panels):
            panel.mu = doublet_strengths[panel_id]
        
        V_fs_norm = self.V_fs.norm()
        
        for panel in mesh.panels:
            
            panel.Velocity = panel_velocity2(panel, mesh.give_neighbours(panel), self.V_fs)            
            
            panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2
    
    def solve2(self, mesh:PanelMesh):
        for panel in mesh.panels:
            panel.sigma = source_strength(panel, self.V_fs) 
        
        B, C = self.influence_coeff_matrices(mesh.panels)
        
        RHS = right_hand_side(mesh.panels, B)
        
        doublet_strengths = np.linalg.solve(C, RHS)
        
        for panel_id, panel in enumerate(mesh.panels):
            panel.mu = doublet_strengths[panel_id]
        
        V_fs_norm = self.V_fs.norm()
        
        for panel in mesh.panels:
            ## 速度计算方法
            # 把所有公共边邻居当中，法向量夹角比较小的三个(相当于是在曲面上还算平滑的几个邻居)取出来
            # 同样是对三个维度进行最小二乘法，但是使用优秀的自己编写的lsq函数。
            # 其实这个方法相对于前面的solve可能相对优胜一点，最起码可以把plane求解的不错。
            # 但是对于其他的问题不好保障欸
            panel_neighbours = mesh.give_neighbours3(panel, 3)
            panel.Velocity = panel_velocity2(panel, panel_neighbours, self.V_fs, rcond = 1e-3)
            
            panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2 
    
    @staticmethod
    def influence_coeff_matrices(panels):
        
        # Katz & Plotkin eq(9.24, 9.25)
        
        n = len(panels)
        B = np.zeros((n, n))
        C = np.zeros_like(B)
        
        for i, panel_i in enumerate(panels):
            
            r_cp = panel_i.r_cp
            
            for j, panel_j in enumerate(panels):
                
                B[i][j], C[i][j] = influence_coeff(r_cp, panel_j)
                
        return B, C
    
    @staticmethod
    def influence_velocoeff_matrices(panels):
        
        # Katz & Plotkin eq(9.24, 9.25)
        
        n = len(panels)
        B = np.zeros((n, n))
        C = np.zeros_like(B)
        
        for i, panel_i in enumerate(panels):
            
            r_cp = panel_i.r_cp
            
            for j, panel_j in enumerate(panels):
                
                B[i][j] = compute_source_panel_velocity(r_cp, panel_j, 1) * panel_i.n
                C[i][j] = compute_dipole_panel_velocity(r_cp, panel_j, 1) * panel_i.n
                
                if i == j:
                    print(B[i][j],C[i][j])
        return B, C


def source_strength(panel, V_fs):
    # Katz & Plotkin eq 9.12
    source_strength = - (panel.n * V_fs)
    return source_strength

def right_hand_side(body_panels, B):

    Nb = len(body_panels)
    RHS = np.zeros(Nb)
    
    for panel_i in body_panels:
        id_i = panel_i.id
        RHS[id_i] = 0
        
        for panel_j in body_panels:
            id_j = panel_j.id
            RHS[id_i] = RHS[id_i] - panel_j.sigma * B[id_i][id_j]
        
    return RHS

def panel_velocity(panel, panel_neighbours, V_fs): 
    # 这里采用的是非精确速度求法，是近似计算，其实可以考虑仿照influence coeff求法求更精确点的(然后炸了)
    # 这个方法其实是忽略了第三个维度，仅仅考虑前两个维度的最小二乘法
    
    n = len(panel_neighbours)
    A = np.zeros((n, 2))
    b = np.zeros((n, 1))
    
    for j, neighbour in enumerate(panel_neighbours):
        r_ij = neighbour.r_cp - panel.r_cp
        r_ij = r_ij.transformation(panel.R)
        A[j][0] = r_ij.x
        A[j][1] = r_ij.y
        b[j][0] = neighbour.mu - panel.mu
    
    del_mu = LeastSquares(A, b)
    components = (-del_mu[0][0], -del_mu[1][0], panel.sigma)
    

    V_disturb = Vector(components)
    V_disturb = V_disturb.transformation(panel.R.T)
    
    V = V_fs + V_disturb
    
    return V

def panel_velocity2(panel, panel_neighbours, V_fs, rcond): 
    # 这里采用的是非精确速度求法，是近似计算，其实可以考虑仿照influence coeff求法求更精确点的
    # 这个方法是考虑了第三个维度的，是三维的最小二乘法。但是用了np内置函数求解lsq
    n = len(panel_neighbours)
    A = np.zeros((n, 3))
    b = np.zeros((n, 1))
    
    for j, neighbour in enumerate(panel_neighbours):
        r_ij = neighbour.r_cp - panel.r_cp
        r_ij = r_ij.transformation(panel.R)
        A[j][0] = r_ij.x
        A[j][1] = r_ij.y
        A[j][2] = r_ij.z
        
        b[j][0] = neighbour.mu - panel.mu
    
    del_mu,_,_,_ = np.linalg.lstsq(A, b, rcond)
    components = (-del_mu[0,0], -del_mu[1,0], panel.sigma)

    V_disturb = Vector(components)
    V_disturb = V_disturb.transformation(panel.R.T)
    
    V = V_fs + V_disturb
    return V

def panel_velocity3(panel, panel_neighbours, V_fs, rcond): 
    # 这里采用的是非精确速度求法，是近似计算，其实可以考虑仿照influence coeff求法求更精确点的
    # 这个方法是考虑了第三个维度的，是三维的最小二乘法。用了自己编写的lsq函数
    n = len(panel_neighbours)
    A = np.zeros((n, 3))
    b = np.zeros((n, 1))
    
    for j, neighbour in enumerate(panel_neighbours):
        r_ij = neighbour.r_cp - panel.r_cp
        r_ij = r_ij.transformation(panel.R)
        A[j][0] = r_ij.x
        A[j][1] = r_ij.y
        A[j][2] = r_ij.z
        
        b[j][0] = neighbour.mu - panel.mu
    
    del_mu = lsq(A, b, rcond)
    components = (-del_mu[0,0], -del_mu[1,0], panel.sigma)

    V_disturb = Vector(components)
    V_disturb = V_disturb.transformation(panel.R.T)
    
    V = V_fs + V_disturb
    return V

def panel_velocity_apame(panel, panel_neighbours, V_fs):
    # 用apame.py里面的函数求解。其实仿照F90写的apame的代码,但是其实效果一般
    ans = cal_velocity(panel, panel_neighbours, 1)
    if ans == False:
        raise Exception("速度计算有错误!")
    a,b = ans
    components = (-a, -b, panel.sigma)

    V_disturb = Vector(components)
    V_disturb = V_disturb.transformation(panel.R.T)
    
    V = V_fs + V_disturb
    return V

def panel_velocity_new(p, mesh, V_fs): 
    # 这里采用的是精确速度求法，但是是错误的！！！
    # 哪里错了呢??
    V_disturb = Vector((0, 0, 0))

    for panel in mesh.panels:
        # 源面板贡献
        V_disturb += compute_source_panel_velocity(p.r_cp, panel, panel.sigma)
        # 双极子面板贡献
        V_disturb += compute_dipole_panel_velocity(p.r_cp, panel, panel.mu)

    return V_fs + V_disturb
   
if __name__ == "__main__":
    from mesh_class import PanelMesh
    from sphere import sphere
    from plot_functions import plot_Cp_SurfaceContours
    
    radius = 1
    num_longitude, num_latitude = 21, 20
    nodes, shells = sphere(radius, num_longitude, num_latitude,
                                     mesh_shell_type='quadrilateral')
    mesh = PanelMesh(nodes, shells)
    V_fs = Vector((1, 0, 0))
    panel_method = Steady_Wakeless_PanelMethod(V_fs)
    panel_method.solve(mesh)
    # print([panel.Cp for panel in mesh.panels])
    plot_Cp_SurfaceContours(mesh.panels)
        