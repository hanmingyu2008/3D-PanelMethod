import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.disturbance_velocity_functions import Dblt_disturb_velocity,Src_disturb_velocity,Vrtx_ring_induced_veloctiy
from PanelMethod3D.influence_coefficient_functions import influence_coeff,Dblt_influence_coeff,compute_source_grouppanel_velocity
from PanelMethod3D.influence_coefficient_functions import compute_source_panel_velocity,compute_dipole_panel_velocity
from PanelMethod3D.mesh_class import PanelMesh,PanelAeroMesh
from PanelMethod3D.LSQ import LeastSquares,lsq
# from PanelMethod3D.apame import cal_velocity
import time
import multiprocessing as mp
from itertools import product
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
        
        B, C = self.influence_coeff_matrices(mesh.panels)

        for panel in mesh.panels:
            panel.sigma = source_strength(panel, self.V_fs) 
        
        RHS = right_hand_side(mesh.panels, B)

        doublet_strengths = np.zeros(mesh.panels_num)
        
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
    
    def solve2(self, mesh:PanelMesh):
        for panel in mesh.panels:
            panel.sigma = source_strength(panel, self.V_fs) 
        time1 = time.time()
        B, C = self.influence_coeff_matrices(mesh.panels)
        time2 = time.time()
        print(f"Constructing Matrix Completed, time cost: {time2 - time1:.2f}",flush=True)
        RHS = right_hand_side(mesh.panels, B)
        print(np.linalg.det(C))
        time1 = time.time()
        doublet_strengths = np.linalg.solve(C, RHS)
        time2 = time.time()
        print(f"Solving Linear Equation Completed, time cost: {time2 - time1:.2f}",flush=True)
        
        for panel_id, panel in enumerate(mesh.panels):
            panel.mu = doublet_strengths[panel_id]
        
        V_fs_norm = self.V_fs.norm()
        time1 = time.time()
        for panel in mesh.panels:
            ## 速度计算方法
            # 把所有公共边邻居当中，法向量夹角比较小的三个(相当于是在曲面上还算平滑的几个邻居)取出来
            # 同样是对三个维度进行最小二乘法，但是使用优秀的自己编写的lsq函数。
            # 其实这个方法相对于前面的solve可能相对优胜一点，最起码可以把plane求解的不错。
            # 但是对于其他的问题不好保障欸
            panel_neighbours = mesh.give_neighbours3(panel, 3)
            panel.Velocity = panel_velocity3(panel, panel_neighbours, self.V_fs, rcond = 1e-2)
            
            panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2 
        time2 = time.time()
        print(f"Computing Velocity Completed, time cost: {time2 - time1:.2f}")

    def solve_newvelo(self, mesh:PanelMesh): # 这个函数还没有改对！！！不要使用！！！
        # 这个函数的B,C和其他的几个并不一样，完全仿照二维的panel method进行的
        # 但是他们不这么干可能也是有他们的道理，这个东西就不指望改对了吧
        # 要说战绩，是能保证panel.Velocity * panel.n = 0,毕竟我方程列的就是这个
        # 但是sphere和解析解对不上啊 
        
        B = self.influence_velocoeff_matrices(mesh.panels)
        
        RHS = np.zeros(len(mesh.panels))

        for panel_id, panel in enumerate(mesh.panels):
            RHS[panel_id] -= panel.n * self.V_fs
        
        source_strengths = np.linalg.solve(B, RHS)
        
        for panel_id, panel in enumerate(mesh.panels):
            panel.sigma = source_strengths[panel_id]
        
        V_fs_norm = self.V_fs.norm()

        for panel in mesh.panels:
            panel.Velocity = panel_velocity_new(panel,mesh,self.V_fs)
            panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm) ** 2 

    def solve_newvelo_group(self, mesh:PanelMesh):

        B,RHS = self.influence_velocoeff_group_matrices_memory(mesh, self.V_fs)
        
        source_strengths = np.linalg.solve(B, RHS)
        
        for panel_id, panel in enumerate(mesh.panels):
            panel.sigma = source_strengths[panel_id]
        
        V_fs_norm = self.V_fs.norm()       

        panel_velocity_group_onlysource(mesh, self.V_fs)
    
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
                
                r_local = (r_cp - panel_j.r_cp).transformation(panel_j.R)
                
        return B, C
    
    @staticmethod
    def influence_velocoeff_matrices(panels):
        
        # Katz & Plotkin eq(9.24, 9.25)
        
        n = len(panels)
        B = np.zeros((n, n))
        
        for i, panel_i in enumerate(panels):
            
            r_cp = panel_i.r_cp
            
            # loop all over panels
            for j, panel_j in enumerate(panels):
                
                B[i][j] = compute_source_panel_velocity(r_cp, panel_j, 1.0) * panel_i.n

        return B
    
    @staticmethod
    def influence_velocoeff_group_matrices_memory(mesh:PanelMesh, V_fs:Vector):
        Coor = np.zeros((3,mesh.panels_num))
        Eta = np.zeros((3,mesh.panels_num))
        for i,panel in enumerate(mesh.panels):
            Coor[0,i] = panel.r_cp.x
            Coor[1,i] = panel.r_cp.y
            Coor[2,i] = panel.r_cp.z
            Eta[0,i] = panel.n.x
            Eta[1,i] = panel.n.y
            Eta[2,i] = panel.n.z
        B = np.zeros((mesh.panels_num,mesh.panels_num))
        for i,panel in enumerate(mesh.panels):
            B[:,i] = np.einsum("ij,ij->j",Eta,compute_source_grouppanel_velocity(Coor,panel))
        
        rhs = np.zeros(mesh.panels_num)
        rhs = - Eta[0] * V_fs.x - Eta[1] * V_fs.y - Eta[2] * V_fs.z
        
        return B,rhs
        
    
class Steady_PanelMethod(PanelMethod):
    
    def __init__(self, V_freestream):
        super().__init__(V_freestream)
    
    def solve(self, mesh:PanelAeroMesh):
            
        body_panels = [mesh.panels[id] for id in mesh.panels_ids["body"]]
        wake_panels = [mesh.panels[id] for id in mesh.panels_ids["wake"]]
        
        for panel in body_panels:
            panel.sigma = source_strength(panel, self.V_fs)
               
        A, B, C = self.influence_coeff_matrices(mesh)
        
        RHS = right_hand_side(body_panels, B)
        
        doublet_strengths = np.linalg.solve(A, RHS)
        
        for panel_i in (body_panels):
            
            panel_i.mu = doublet_strengths[panel_i.id]
            
            if panel_i.id in mesh.TrailingEdge["suction side"]:
                for id_j in mesh.wake_sheddingPanels[panel_i.id]:
                    panel_j = mesh.panels[id_j]
                    panel_j.mu = panel_j.mu + doublet_strengths[panel_i.id]
                    
            elif panel_i.id in mesh.TrailingEdge["pressure side"]:
                for id_j in mesh.wake_sheddingPanels[panel_i.id]:
                    panel_j = mesh.panels[id_j]
                    panel_j.mu = panel_j.mu - doublet_strengths[panel_i.id]
        
        V_fs_norm = self.V_fs.norm()

        ## 一下是一个希腊老哥写的代码注释，其中希腊语汉译如下:
        
        # 这个方法不起作用。我不明白为什么
        # 它更直接（与上述需要近似求解强度μ的梯度的方法相反）
        # 尽管它要慢得多
        
        # for panel in body_panels:
            
        #     # Velocity caclulation with least squares approach (faster)
            
        #     panel_neighbours = mesh.give_neighbours(panel)
        #     panel.Velocity = panel_velocity(panel, panel_neighbours, self.V_fs)
            
          
        #     # Velocity calculation with disturbance velocity functions
            
        #     # Δεν δουλεύει αυτή η μέθοδος. Δεν μπορώ να καταλάβω γιατ΄ί
        #     # Είναι πιο straight forward (σε αντίθεση με την παραπάνω μέθοδο
        #     # που απαιτεί προσεγγιστική επίλυση των gradients της έντασης μ)
        #     # παρ' ότι πολύ πιο αργή
            
        #     # panel.Velocity = Velocity(mesh.V_fs, panel.r_cp, body_panels,
        #     #                           wake_panels)
            
            
        #     # pressure coefficient calculation
        #     panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2

        # 还是采用基于 "Program VSAERO Theory Document" 的方法为好,似乎也不怎么样
        # VSAERO_onbody_analysis(self.V_fs, mesh)
        for panel_id in mesh.panels_ids["body"]:
            panel = mesh.panels[panel_id]
            '''
            panel_neighbours = mesh.give_neighbours3(panel, 3)
            panel.Velocity = panel_velocity2(panel, panel_neighbours, self.V_fs, rcond = 1e-3)
            '''
            panel.Velocity = panel_velocity(panel,mesh.give_neighbours(panel),self.V_fs)
            panel.Cp = 1 - (panel.Velocity.norm()/V_fs_norm)**2
        
    @staticmethod
    def influence_coeff_matrices(mesh:PanelAeroMesh):
        
        Nb = len(mesh.panels_ids["body"])
        Nw = len(mesh.panels_ids["wake"])
        B = np.zeros((Nb, Nb))
        C = np.zeros((Nb, Nb+Nw))
        A = np.zeros_like(B)
        
        for id_i in mesh.panels_ids["body"]:
            panel_i = mesh.panels[id_i]
            r_cp = panel_i.r_cp
            
            for id_j in mesh.panels_ids["body"]:
                
                panel_j = mesh.panels[id_j]

                B[id_i][id_j], C[id_i][id_j] = influence_coeff(r_cp, panel_j)
                A[id_i][id_j] = C[id_i][id_j]
                
                if id_j in mesh.TrailingEdge["suction side"]:
                    for id_k in mesh.wake_sheddingPanels[id_j]:
                        panel_k = mesh.panels[id_k]
                        C[id_i][id_k] = Dblt_influence_coeff(r_cp, panel_k)
                        A[id_i][id_j] = A[id_i][id_j] + C[id_i][id_k]
                        
                elif id_j in mesh.TrailingEdge["pressure side"]:
                    for id_k in mesh.wake_sheddingPanels[id_j]:
                        panel_k = mesh.panels[id_k]
                        C[id_i][id_k] = Dblt_influence_coeff(r_cp, panel_k)
                        A[id_i][id_j] = A[id_i][id_j] - C[id_i][id_k]
                
        return A, B, C


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

def panel_velocity_new(p, mesh, V_fs): 
    # 这里采用的是精确速度求法，但是是错误的！！！
    # 哪里错了呢??
    V_disturb = Vector((0, 0, 0))

    for panel in mesh.panels:
        # 源面板贡献
        if panel.sigma != 0:
            V_disturb += compute_source_panel_velocity(p.r_cp, panel, panel.sigma)
        # 双极子面板贡献
        if panel.mu != 0:
            V_disturb += compute_dipole_panel_velocity(p.r_cp, panel, panel.mu)

    return V_fs + V_disturb


def panel_velocity_group_onlysource(mesh:PanelMesh, V_fs): 
    # 这里采用的是精确速度求法，但是是错误的！！！
    # 哪里错了呢??
    V_disturb = np.zeros((3,mesh.panels_num))
    V_disturb[0] += V_fs.x
    V_disturb[1] += V_fs.y
    V_disturb[2] += V_fs.z

    Coor = np.zeros_like(V_disturb)
    for i,panel in enumerate(mesh.panels):
            Coor[0,i] = panel.r_cp.x
            Coor[1,i] = panel.r_cp.y
            Coor[2,i] = panel.r_cp.z

    for panel in mesh.panels:
            
            V_disturb += panel.sigma * compute_source_grouppanel_velocity(Coor,panel)

    for i in range(mesh.panels_num):
        mesh.panels[i].Velocity = Vector((V_disturb[0,i],V_disturb[1,i],V_disturb[2,i]))    
        mesh.panels[i].Cp = 1 - (mesh.panels[i].Velocity.norm()/V_fs.norm())**2   