import numpy as np
from vector_class import Vector
from panel_class import Panel, quadPanel, triPanel
from is_inside_polygon import is_inside_polygon


def Src_influence_coeff(r_p:Vector, panel:Panel, alpha=10):
    n = panel.num_vertices
    r_vertex = panel.r_vertex_local
    r_cp = panel.r_cp
    R = panel.R 
    
    r = r_p - r_cp
    r_local = r.transformation(R)
    r = r_local
    
    if r.norm() >= alpha * panel.char_length:
        
        B = panel.area/r.norm()
    
    elif r.z == 0:
        
        B = 0
        for i in range(n-1, -1, -1):
            # panel numbering follow counter clock wise direction
            # Hess and Smith integrals are calculated with clock wise ordering 
            a = (i+1)%n  # 0, 3, 2, 1 (cw) (instead 0, 1, 2, 3 (cw))
            b = i # 3, 2, 1, 0, (clock wise) (instead 1, 2, 3, 0 (cw))
            
            r_ab = r_vertex[b] - r_vertex[a]
            d_ab = r_ab.norm()
            r_a = r - r_vertex[a]
            r_a = r_a.norm()
            r_b = r - r_vertex[b]
            r_b = r_b.norm()     
            
            first_term = (
                (r.x - r_vertex[a].x) * (r_vertex[b].y - r_vertex[a].y)
                - (r.y - r_vertex[a].y) * (r_vertex[b].x - r_vertex[a].x)
                )
            
            first_term = first_term/d_ab
               
            if (r_a + r_b - d_ab) == 0:
                # point p coincide lies on a panel's edge
                # log term ---> inf => B ---> inf
                log_term = 0
        
            else:
                # katz & Plotkin
                log_term = np.log((r_a + r_b + d_ab)/(r_a + r_b - d_ab))
                
                # paper of Lothar birk  
                # log_term = np.log((r_a + r_b - d_ab)/(r_a + r_b + d_ab))
                
            B = B + first_term * log_term  
              
    else:
        
        B = 0
        for i in range(n-1, -1, -1):
            # panel numbering follow counter clock wise direction
            # Hess and Smith integrals are calculated with clock wise ordering 
            a = (i+1)%n  # 0, 3, 2, 1 (cw) (instead 0, 1, 2, 3 (cw))
            b = i # 3, 2, 1, 0, (clock wise) (instead 1, 2, 3, 0 (cw))
            
            r_ab = r_vertex[b] - r_vertex[a]
            d_ab = r_ab.norm()
            r_a = r - r_vertex[a]
            r_a = r_a.norm()
            r_b = r - r_vertex[b]
            r_b = r_b.norm()       
            
            first_term = (
                (r.x - r_vertex[a].x) * (r_vertex[b].y - r_vertex[a].y)
                - (r.y - r_vertex[a].y) * (r_vertex[b].x - r_vertex[a].x)
                )
            
            first_term = first_term/d_ab
            
            # katz & Plotkin
            log_term = np.log((r_a + r_b + d_ab)/(r_a + r_b - d_ab))
            
            # paper of Lothar birk  
            # log_term = np.log((r_a + r_b - d_ab)/(r_a + r_b + d_ab))
                
            if (r_vertex[b].x - r_vertex[a].x) == 0:
                m_ab = np.inf 
            else:
                m_ab = ( (r_vertex[b].y - r_vertex[a].y)
                        /(r_vertex[b].x - r_vertex[a].x) )
            
            e_a = (r.x - r_vertex[a].x)**2 + r.z**2
            e_b = (r.x - r_vertex[b].x)**2 + r.z**2
            h_a = (r.x - r_vertex[a].x)*(r.y - r_vertex[a].y)
            h_b = (r.x - r_vertex[b].x)*(r.y - r_vertex[b].y)
            
            arctan_term = ( np.arctan((m_ab*e_a - h_a)/(r.z*r_a))
                           - np.arctan((m_ab*e_b - h_b)/(r.z*r_b)) )
            
            # Katz & Plotkin            
            # B = B + first_term * log_term - np.abs(r.z) * arctan_term
            
            # paper of Lothar birk            
            B = B + first_term * log_term - r.z * arctan_term
                  
    B = - 1/(4 * np.pi) * B
    
    return B

def Dblt_influence_coeff(r_p:Vector, panel:Panel, alpha=10):
    n = panel.num_vertices
    r_vertex = panel.r_vertex_local
    r_cp = panel.r_cp
    R = panel.R 
    
    r = r_p - r_cp
    r_local = r.transformation(R)
    r = r_local
    
    if r.norm() >= alpha * panel.char_length:
        C = - 1/(4*np.pi) * ( panel.area * r.z)/(r.norm()**3)
    
    elif r.z == 0:
        point = (r.x, r.y)  
        polygon = [(r_vertex[i].x, r_vertex[i].y) for i in range(n)]
        
        if is_inside_polygon(polygon, point):
            # point p lies on panel's surface as z-->0
            # C = - 0.5  # if z--> +0
            C = - 0.5  # if z--> -0
        else:
            # point p lies outside of panel's surface as z-->0
            C = 0
        
    else:
        C = 0
        for i in range(n-1, -1, -1):
            # panel numbering follow counter clock wise direction
            # Hess and Smith integrals are calculated with clock wise ordering 
            a = (i+1)%n  # 0, 3, 2, 1 (cw) (instead 0, 1, 2, 3 (cw))
            b = i # 3, 2, 1, 0, (clock wise) (instead 1, 2, 3, 0 (cw))
            
            r_a = r - r_vertex[a]
            r_a = r_a.norm()
            r_b = r - r_vertex[b]
            r_b = r_b.norm()   
            
            if (r_vertex[b].x - r_vertex[a].x) == 0:
                m_ab = np.inf 
            else:
                m_ab = ( (r_vertex[b].y - r_vertex[a].y)
                        /(r_vertex[b].x - r_vertex[a].x) )
            
            e_a = (r.x - r_vertex[a].x)**2 + r.z**2
            e_b = (r.x - r_vertex[b].x)**2 + r.z**2
            h_a = (r.x - r_vertex[a].x)*(r.y - r_vertex[a].y)
            h_b = (r.x - r_vertex[b].x)*(r.y - r_vertex[b].y)
            
            arctan_term = ( np.arctan((m_ab*e_a - h_a)/(r.z*r_a))
                           - np.arctan((m_ab*e_b - h_b)/(r.z*r_b)) )
                
            C = C + arctan_term
            
        C = - 1/(4*np.pi) * C
        
    return C

def influence_coeff(r_p:Vector, panel:Panel, alpha=10):

    ## 参见 Lothar Birk 的文章。

    n = panel.num_vertices
    r_vertex = panel.r_vertex_local
    r_cp = panel.r_cp
    R = panel.R 
    
    r = r_p - r_cp
    r_local = r.transformation(R)
    r = r_local
    
    # if r.norm() >= alpha * panel.char_length:
    #    B = panel.area/r.norm()
    #    C = ( panel.area * r.z)/(r.norm()**3)
    
    if r.z == 0:
        
        point = (r.x, r.y)  
        polygon = [(r_vertex[i].x, r_vertex[i].y) for i in range(n)]
    
        if is_inside_polygon(polygon, point):
            # C =  2 * np.pi  # if z--> +0
            C = -2 * np.pi  # if z--> -0
        else:
            C = 0
        
        B = 0
        for i in range(n-1, -1, -1):
            # 我这个panel它是逆时针。
            # Hess-Smith 积分是顺时针。 
            a = (i+1)%n  # 0, 3, 2, 1 (cw) (instead 0, 1, 2, 3 (cw))
            b = i # 3, 2, 1, 0, (clock wise) (instead 1, 2, 3, 0 (cw))
            
            r_ab = r_vertex[b] - r_vertex[a]
            d_ab = r_ab.norm()
            r_a = r - r_vertex[a]
            r_a = r_a.norm()
            r_b = r - r_vertex[b]
            r_b = r_b.norm()     
            
            first_term = (
                (r.x - r_vertex[a].x) * (r_vertex[b].y - r_vertex[a].y)
                - (r.y - r_vertex[a].y) * (r_vertex[b].x - r_vertex[a].x)
                )
            
            first_term = first_term/d_ab
            
            if (r_a + r_b - d_ab) == 0:
                # point p coincide lies on a panel's edge
                # log term ---> inf => B ---> inf
                log_term = 0
        
            else:
                # katz & Plotkin
                log_term = np.log((r_a + r_b + d_ab)/(r_a + r_b - d_ab))
                
                # paper of Lothar birk  
                # log_term = np.log((r_a + r_b - d_ab)/(r_a + r_b + d_ab))
                
            B = B + first_term * log_term  
                           
    else:
        B = 0
        C = 0
        for i in range(n-1, -1, -1):
            # panel numbering follow counter clock wise direction
            # Hess and Smith integrals are calculated with clock wise ordering 
            a = (i+1)%n  # 0, 3, 2, 1 (cw) (instead 0, 1, 2, 3 (cw))
            b = i # 3, 2, 1, 0, (clock wise) (instead 1, 2, 3, 0 (cw))
            
            r_ab = r_vertex[b] - r_vertex[a]
            d_ab = r_ab.norm()
            r_a = r - r_vertex[a]
            r_a = r_a.norm()
            r_b = r - r_vertex[b]
            r_b = r_b.norm()          
            
            first_term = (
                (r.x - r_vertex[a].x) * (r_vertex[b].y - r_vertex[a].y)
                - (r.y - r_vertex[a].y) * (r_vertex[b].x - r_vertex[a].x)
                )
            
            first_term = first_term/d_ab
            
            # katz & Plotkin
            log_term = np.log((r_a + r_b + d_ab)/(r_a + r_b - d_ab))
            
            if (r_vertex[b].x - r_vertex[a].x) == 0:
                m_ab = np.inf 
            else:
                m_ab = ( (r_vertex[b].y - r_vertex[a].y)
                        /(r_vertex[b].x - r_vertex[a].x) )
            
            e_a = (r.x - r_vertex[a].x)**2 + r.z**2
            e_b = (r.x - r_vertex[b].x)**2 + r.z**2
            h_a = (r.x - r_vertex[a].x)*(r.y - r_vertex[a].y)
            h_b = (r.x - r_vertex[b].x)*(r.y - r_vertex[b].y)
            
             
            arctan_term = ( np.arctan((m_ab*e_a - h_a)/(r.z*r_a))
                           - np.arctan((m_ab*e_b - h_b)/(r.z*r_b)) )
     
            
            
            # Katz & Plotkin            
            # B = B + first_term * log_term - np.abs(r.z) * arctan_term
            
            # paper of Lothar birk            
            B = B + first_term * log_term - r.z * arctan_term
            
            # Katz & Plotkin  
            C = C + arctan_term 
            
    B = - 1/(4 * np.pi) * B
    C = - 1/(4 * np.pi) * C
        
    return B, C

import numpy as np
from math import sqrt, log, atan, tan

def compute_source_panel_velocity(p_g: Vector, panel: Panel, sigma: float) -> Vector:
    """
    计算恒定源强度面板对场点p的速度影响
    
    参数:
        p_g: 场点全局坐标 (Vector)
        panel: 面板对象 (Panel或其子类)
        sigma: 面板源强度
        
    返回:
        v_g: 全局坐标系中的速度向量 (Vector)
    """
    # 1. 将场点转换到面板局部坐标系
    p_local = (p_g - panel.r_cp).transformation(panel.R)
    x, y, z = p_local.x, p_local.y, p_local.z

    if abs(z) < 1e-10:
        return 0.5 * sigma * panel.n
    
    # 2. 初始化速度分量
    phi_x = 0.0
    phi_y = 0.0
    phi_z = 0.0
    
    # 3. 对每条边计算贡献 (论文中的公式47-49)
    for k in range(panel.num_vertices):
        # 当前边端点索引
        k_next = (k - 1) % panel.num_vertices
        
        # 获取局部坐标中的顶点
        q_k = panel.r_vertex_local[k]
        q_k_next = panel.r_vertex_local[k_next]
        
        x_k, y_k = q_k.x, q_k.y
        x_k_next, y_k_next = q_k_next.x, q_k_next.y
        
        # 计算边长度
        d_k = sqrt((x_k_next - x_k)**2 + (y_k_next - y_k)**2)
        if d_k == 0:  # 跳过零长度边
            continue
        
        # 计算距离r_k和r_k+1
        r_k = sqrt((x - x_k)**2 + (y - y_k)**2 + z**2)
        r_k_next = sqrt((x - x_k_next)**2 + (y - y_k_next)**2 + z**2)
        
        # 计算对数项
        log_term = log((r_k + r_k_next - d_k) / (r_k + r_k_next + d_k))
        
        # x方向分量 (公式47)
        phi_x += (y_k_next - y_k) * log_term / d_k 
        
        # y方向分量 (公式48)
        phi_y += -(x_k_next - x_k) * log_term / d_k 
        
        # z方向分量 (公式49)
        m_k = (y_k_next - y_k) / (x_k_next - x_k) if (x_k_next - x_k) != 0 else float('inf')
        e_k = (x - x_k)**2 + z**2
        h_k = (x - x_k) * (y - y_k)
        e_k_next = (x - x_k_next)**2 + z**2
        h_k_next = (x - x_k_next) * (y - y_k_next)
        
        atan_term1 = atan((m_k * e_k - h_k) / (z * r_k)) if z != 0 else 0
        atan_term2 = atan((m_k * e_k_next - h_k_next) / (z * r_k_next)) if z != 0 else 0
        phi_z += atan_term1 - atan_term2
    
    # 4. 乘以系数和源强度
    factor = -sigma / (4 * np.pi)
    phi_x *= factor
    phi_y *= factor
    phi_z *= factor
    
    # 5. 将速度转换回全局坐标系
    v_local = Vector((phi_x, phi_y, phi_z))
    v_g = v_local.transformation(panel.R.T)  # 转置矩阵用于反向变换
    
    return v_g

def compute_dipole_panel_velocity(p_g: Vector, panel: Panel, mu: float) -> Vector:
    """
    计算恒定双极子强度面板对场点p的速度影响
    
    参数:
        p_g: 场点全局坐标 (Vector)
        panel: 面板对象 (Panel或其子类)
        mu: 面板双极子强度
        
    返回:
        v_g: 全局坐标系中的速度向量 (Vector)
    """
    # 1. 将场点转换到面板局部坐标系
    p_local = (p_g - panel.r_cp).transformation(panel.R)
    x, y, z = p_local.x, p_local.y, p_local.z

    if abs(z) < 1e-10:
        return Vector((0,0,0))
    
    # 2. 初始化速度分量
    psi_x = 0.0
    psi_y = 0.0
    psi_z = 0.0
    
    # 3. 对每条边计算贡献 (论文中的公式60和52,55,58)
    for k in range(panel.num_vertices):
        # 当前边端点索引
        k_next = (k - 1) % panel.num_vertices
        
        # 获取局部坐标中的顶点
        q_k = panel.r_vertex_local[k]
        q_k_next = panel.r_vertex_local[k_next]
        
        x_k, y_k = q_k.x, q_k.y
        x_k_next, y_k_next = q_k_next.x, q_k_next.y
        
        # 计算边长度
        d_k = sqrt((x_k_next - x_k)**2 + (y_k_next - y_k)**2)
        if d_k == 0:  # 跳过零长度边
            continue
        
        # 计算距离r_k和r_k+1
        r_k = sqrt((x - x_k)**2 + (y - y_k)**2 + z**2)
        r_k_next = sqrt((x - x_k_next)**2 + (y - y_k_next)**2 + z**2)
        
        # 计算辅助变量 (公式43-44)
        rho_k = r_k * r_k_next + (x - x_k)*(x - x_k_next) + (y - y_k)*(y - y_k_next) + z**2
        lambda_k = (x - x_k)*(y - y_k_next) - (x - x_k_next)*(y - y_k)
        
        # x方向分量 (公式52的负值)
        term_x = z * (y_k_next - y_k) * (r_k + r_k_next) / (r_k * r_k_next * rho_k)
        psi_x -= term_x
        
        # y方向分量 (公式55的负值)
        term_y = -z * (x_k_next - x_k) * (r_k + r_k_next) / (r_k * r_k_next * rho_k)
        psi_y -= term_y
        
        # z方向分量 (公式58的负值)
        term_z = lambda_k * (r_k + r_k_next) / (r_k * r_k_next * rho_k)
        psi_z -= term_z
    
    # 4. 乘以系数和双极子强度
    factor = -mu / (4 * np.pi)
    psi_x *= factor
    psi_y *= factor
    psi_z *= factor
    
    # 5. 将速度转换回全局坐标系
    v_local = Vector((psi_x, psi_y, psi_z))
    v_g = v_local.transformation(panel.R.T)  # 转置矩阵用于反向变换
    
    return v_g

if __name__=='__main__':
    from matplotlib import pyplot as plt
    vertex1 = Vector((-1, -1, 1))
    vertex2 = Vector((1, -1, 1))
    vertex3 = Vector((1, 1, 1))
    vertex4 = Vector((-1, 1, 1))
    
    # Quadrilateral panel
    panel = quadPanel(vertex1, vertex2, vertex3, vertex4)
    
    # Triangular panel
    # panel = triPanel(vertex1, vertex2, vertex3)
    
    r_p = panel.r_cp

    C = Dblt_influence_coeff(r_p, panel)
    B = Src_influence_coeff(r_p, panel)
    print("C = " + str(C) + "\nB = " + str(B))
    
    B, C = influence_coeff(r_p, panel)
    print("C = " + str(C) + "\nB = " + str(B))
    
    
    # plot panel
    
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
  