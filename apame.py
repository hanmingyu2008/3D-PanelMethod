import numpy as np
from vector_class import Vector
from mesh_class import PanelMesh
from scipy.linalg import solve, lstsq

def cal_velocity(panel, neighboor_list, order):

    coplanarity_angle = 0.15
    if len(neighboor_list) == 4:
        if order == 2:
            ql, qp, v_err = quadratic_4(panel, neighboor_list, coplanarity_angle)
        else:
            ql, qp, v_err = linear_4(panel, neighboor_list, coplanarity_angle)
    elif len(neighboor_list) == 3:
        if order == 2:
            ql, qp, v_err = quadratic_3(panel, neighboor_list, coplanarity_angle)
        else:
            ql, qp, v_err = linear_3(panel, neighboor_list, coplanarity_angle)
    else:
        return False
    if v_err :
        return False
    return ql, qp
    

    # ====================== 插值核心函数 ======================

def quadratic_4(panel, neighboor_list, coplanarity_angle):
        """
        Calculate induced velocities using quadratic interpolation with 5 points
        """
        T = np.zeros((4, 4))
        mu = np.zeros((4,1))
        
        for i,neigh in enumerate(neighboor_list):
            radius_vector = neigh.r_cp-panel.r_cp
            radius_vector = radius_vector.transformation(panel.R)
            correcting_factor = curvature_correction(coplanarity_angle, radius_vector, panel.n, neigh.n)
            
            T[i, 0] = radius_vector.x * correcting_factor
            T[i, 1] = radius_vector.y * correcting_factor
            T[i, 2] = T[i, 0]**2
            T[i, 3] = T[i, 1]**2

            mu[i] = neigh.mu - panel.mu
        
        try:
            solution, _, _, _ = lstsq(T, mu, cond=1e-5)
            ql = solution[0,0]  # A coefficient
            qp = solution[1,0]  # B coefficient
            v_err = False
        except:
            ql = 0
            qp = 0
            v_err = True

        return ql,qp,v_err
    
def linear_4(panel, neighboor_list, coplanarity_angle):
        """
        Calculate induced velocities using linear interpolation with 5 points
        """
        T = np.zeros((4, 2))
        mu = np.zeros((4,1))
        
        for i,neigh in enumerate(neighboor_list):
            radius_vector = neigh.r_cp-panel.r_cp
            radius_vector = radius_vector.transformation(panel.R)
            correcting_factor = curvature_correction(coplanarity_angle, radius_vector, panel.n, neigh.n)
            
            T[i, 0] = radius_vector.x * correcting_factor
            T[i, 1] = radius_vector.y * correcting_factor

            mu[i] = neigh.mu - panel.mu
        
        try:
            solution, _, _, _ = lstsq(T, mu, cond=1e-5)
            ql = solution[0,0]  # A coefficient
            qp = solution[1,0]  # B coefficient
            v_err = False
        except:
            ql = 0
            qp = 0
            v_err = True

        return ql,qp,v_err
    
def quadratic_3(panel, neighboor_list, coplanarity_angle):
        """
        Calculate induced velocities using quadratic interpolation with 4 points
        """
        T = np.zeros((3, 3))
        mu = np.zeros((3,1))
        
        for i,neigh in enumerate(neighboor_list):
            radius_vector = neigh.r_cp-panel.r_cp
            radius_vector = radius_vector.transformation(panel.R)
            correcting_factor = curvature_correction(coplanarity_angle, radius_vector, panel.n, neigh.n)
            
            T[i, 0] = radius_vector.x * correcting_factor
            T[i, 1] = radius_vector.y * correcting_factor
            T[i, 2] = T[i, 0]**2+T[i, 1]**2

            mu[i] = neigh.mu - panel.mu
        
        try:
            solution, _, _, _ = lstsq(T, mu, cond=1e-5)
            ql = solution[0,0]  # A coefficient
            qp = solution[1,0]  # B coefficient
            v_err = False
        except:
            ql = 0
            qp = 0
            v_err = True

        return ql,qp,v_err

def linear_3(panel, neighboor_list, coplanarity_angle):
        """
        四点线性插值
        参数：
            points: 四个点的坐标 [(x1,y1,z1), ...]
            l_vec: 面板纵向单位向量
            p_vec: 面板横向单位向量
            normals: 四个点的法向量
            gama_nodal: 节点双偶极子强度
        返回：
            ql, qp: 局部坐标系速度分量
            v_err: 错误标志
        """
        T = np.zeros((3, 2))
        mu = np.zeros((3,1))
        
        for i,neigh in enumerate(neighboor_list):
            radius_vector = neigh.r_cp-panel.r_cp
            radius_vector = radius_vector.transformation(panel.R)
            correcting_factor = curvature_correction(coplanarity_angle, radius_vector, panel.n, neigh.n)
            
            T[i, 0] = radius_vector.x * correcting_factor
            T[i, 1] = radius_vector.y * correcting_factor

            mu[i] = neigh.mu - panel.mu                    

        try:
            solution, _, _, _ = lstsq(T, mu, cond=1e-5)
            ql = solution[0,0]  # A coefficient
            qp = solution[1,0]  # B coefficient
            v_err = False
        except:
            ql = 0
            qp = 0
            v_err = True
        
        return ql,qp,v_err

def curvature_correction(coplanarity_angle, radius_vector, normal_1, normal_2):
        """
        曲率修正因子计算
        参数：
            r_vec: 半径向量 (从当前面板指向相邻面板)
            n1, n2: 两个面板的法向量
            angle_thresh: 共面性角度阈值(弧度)
        返回：
            修正因子 (float)
        """
        pi = 2.0 * np.arccos(0.0)
        
        # Calculate radius vector length
        side_c = np.linalg.norm(radius_vector)
        radius_unit_vector = radius_vector / side_c
        
        # Calculate average normal
        average_normal = (normal_1 + normal_2) / 2
        average_normal_mod = np.linalg.norm(average_normal)
        average_normal = average_normal / average_normal_mod
        
        # Calculate projection plane
        plane_unit_vector = Vector.cross_product(radius_unit_vector, average_normal)
        
        # Project normals to projection plane
        scalar = normal_2 * plane_unit_vector
        normal_2_proj = normal_2 - scalar * plane_unit_vector
        scalar = normal_1 * plane_unit_vector
        normal_1_proj = normal_1 - scalar * plane_unit_vector
        
        # Calculate angles of triangle
        temp = (radius_unit_vector * normal_2_proj) / (np.linalg.norm(radius_unit_vector) * np.linalg.norm(normal_2_proj))
        if temp > 1: temp = 1
        if temp < -1: temp = -1
        ang1 = abs(np.arccos(temp) - pi/2)
        temp = (radius_unit_vector * normal_1_proj) / (np.linalg.norm(radius_unit_vector) * np.linalg.norm(normal_1_proj))
        if temp > 1: temp = 1
        if temp < -1: temp = -1
        ang2 = abs(np.arccos(temp) - pi/2)
        ang3 = pi - ang1 - ang2
        
        if ang3 < coplanarity_angle and ang3 > pi - coplanarity_angle:
            # Calculate triangle sides using sine law
            return abs((np.sin(ang1) + np.sin(ang2)) / np.sin(ang3))
        else:
            return 1.0

    # ====================== 辅助函数 ======================

def handle_error(self, method, speed):
        """处理插值错误"""
        if method == 0:
            return np.zeros_like(speed), np.zeros_like(speed)
        else:
            return np.sqrt(speed), np.sqrt(speed)