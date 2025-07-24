import numpy as np
from vector_class import Vector
from panel_class import Panel, quadPanel, triPanel
from is_inside_polygon import is_inside_polygon
# 几个函数都是参见Lothar的论文,第一个是对速度势的贡献,后两个是对速度的贡献。
# 这里“贡献”指的都是单位大小的source/doublet面板对某个点处速度(势)的影响。

def Dblt_influence_coeff(p_g:Vector, panel:Panel):
    
    C = 0
    for k in range(panel.num_vertices):
        k_next = (k+1)%(panel.num_vertices)
        q_k = panel.r_vertex_local[k]
        q_k_next = panel.r_vertex_local[k_next]

        p_local = (p_g - panel.r_cp).transformation(panel.R)
        x, y, z = p_local.x, p_local.y, p_local.z
        
        x_k, y_k = q_k.x, q_k.y
        x_k_next, y_k_next = q_k_next.x, q_k_next.y

        d_k = sqrt((x_k_next-x_k)**2+(y_k_next-y_k)**2)
        if d_k == 0:
            continue
        e_k = (x-x_k)**2 + z**2
        e_k_next = (x-x_k_next)**2 + z**2
        h_k = (x-x_k)*(y-y_k)
        h_k_next = (x-x_k_next)*(y-y_k_next)
        m_k = (y_k_next-y_k)/(x_k_next-x_k) if x_k != x_k_next else float("inf")
        r_k = sqrt((x-x_k)**2+(y-y_k)**2+z**2)
        r_k_next = sqrt((x-x_k_next)**2+(y-y_k_next)**2+z**2)

        termc = np.arctan((m_k*e_k-h_k)/(z*r_k)) - np.arctan((m_k*e_k_next-h_k_next)/(z*r_k_next)) if z!=0 else 0

        C -= termc 
    
    if z == 0:
        
        point = (x,y)  
        polygon = [(panel.r_vertex_local[i].x, panel.r_vertex_local[i].y) for i in range(panel.num_vertices)]
    
        if is_inside_polygon(polygon, point):
            # C =  2 * np.pi  # if z--> +0
            C = -2 * np.pi  # if z--> -0
        else:
            C = 0
    
    C = - 1/(4 * np.pi) * C
        
    return C

def influence_coeff(p_g:Vector, panel:Panel):
    B = 0
    C = 0
    for k in range(panel.num_vertices):
        k_next = (k+1)%(panel.num_vertices)
        q_k = panel.r_vertex_local[k]
        q_k_next = panel.r_vertex_local[k_next]

        p_local = (p_g - panel.r_cp).transformation(panel.R)
        x, y, z = p_local.x, p_local.y, p_local.z
        
        x_k, y_k = q_k.x, q_k.y
        x_k_next, y_k_next = q_k_next.x, q_k_next.y

        d_k = sqrt((x_k_next-x_k)**2+(y_k_next-y_k)**2)
        if d_k == 0:
            continue
        e_k = (x-x_k)**2 + z**2
        e_k_next = (x-x_k_next)**2 + z**2
        h_k = (x-x_k)*(y-y_k)
        h_k_next = (x-x_k_next)*(y-y_k_next)
        m_k = (y_k_next-y_k)/(x_k_next-x_k) if x_k != x_k_next else float("inf")
        r_k = sqrt((x-x_k)**2+(y-y_k)**2+z**2)
        r_k_next = sqrt((x-x_k_next)**2+(y-y_k_next)**2+z**2)

        termb1 = ((x-x_k)*(y_k_next-y_k)-(y-y_k)*(x_k_next-x_k)) * np.log((r_k+r_k_next-d_k)/(r_k+r_k_next+d_k)) / d_k
        termb2 = z*(np.arctan((m_k*e_k-h_k)/(z*r_k)) - np.arctan((m_k*e_k_next-h_k_next)/(z*r_k_next))) if z!=0 else 0
        termc = np.arctan((m_k*e_k-h_k)/(z*r_k)) - np.arctan((m_k*e_k_next-h_k_next)/(z*r_k_next)) if z!=0 else 0

        B += termb1 + termb2
        C -= termc 
    
    if z == 0:
        
        point = (x,y)  
        polygon = [(panel.r_vertex_local[i].x, panel.r_vertex_local[i].y) for i in range(panel.num_vertices)]
    
        if is_inside_polygon(polygon, point):
            # C =  2 * np.pi  # if z--> +0
            C = -2 * np.pi  # if z--> -0
        else:
            C = 0
    
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
    '''
    if z == 0:
        point = (x,y)  
        polygon = [(panel.r_vertex_local[i].x, panel.r_vertex_local[i].y) for i in range(panel.num_vertices)]
    
        if is_inside_polygon(polygon, point):
            return 0.5 * sigma * panel.n
        else:
            return Vector((0,0,0))'''
    
    # 2. 初始化速度分量
    phi_x = 0.0
    phi_y = 0.0
    phi_z = 0.0
    
    # 3. 对每条边计算贡献 (论文中的公式47-49)
    for k in range(panel.num_vertices):
        # 当前边端点索引
        k_next = (k + 1) % panel.num_vertices
        
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
        log_term = np.log((r_k + r_k_next - d_k) / (r_k + r_k_next + d_k))
        
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
        
        atan_term1 = np.arctan((m_k * e_k - h_k) / (z * r_k)) if z != 0 else np.pi/2
        atan_term2 = np.arctan((m_k * e_k_next - h_k_next) / (z * r_k_next)) if z != 0 else np.pi/2
        phi_z += atan_term1 - atan_term2

    if z == 0 :
        point = (x,y)  
        polygon = [(panel.r_vertex_local[i].x, panel.r_vertex_local[i].y) for i in range(panel.num_vertices)]
    
        if is_inside_polygon(polygon, point):
            phi_z = -2*np.pi
        else:
            phi_z = 0

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
    
    # 2. 初始化速度分量
    psi_x = 0.0
    psi_y = 0.0
    psi_z = 0.0
    
    # 3. 对每条边计算贡献 (论文中的公式60和52,55,58)
    for k in range(panel.num_vertices):
        # 当前边端点索引
        k_next = (k + 1) % panel.num_vertices
        
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
  