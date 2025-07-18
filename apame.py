import numpy as np
from vector_class import Vector
from mesh_class import PanelMesh
from scipy.linalg import solve, lstsq

def cal_velocity(panel, neighboor_list, order):
    if len(neighboor_list) == 4:
        if order == 2:
            ql, qp, v_err = quadratic_4(panel, neighboor_list)
        else:
            ql, qp, v_err = linear_4(panel, neighboor_list)
    elif len(neighboor_list) == 3:
        if order == 2:
            ql, qp, v_err = quadratic_3(panel, neighboor_list)
        else:
            ql, qp, v_err = linear_3(panel, neighboor_list)
    else:
        return False
    if v_err :
        return False
    

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
            T[i, 4] = 1
            mu[i] = neigh.mu - panel.mu
        
        try:
            solution = solve(T, mu)
            ql = solution[0]  # A coefficient
            qp = solution[1]  # B coefficient
            v_err = False
        except:
            ql = 0
            qp = 0
            v_err = 1

        return ql,qp,v_err
    
def linear_4(self, x, y, z, l, p, n, coplanarity_angle, gama):
        """
        Calculate induced velocities using linear interpolation with 5 points
        """
        T = np.zeros((5, 3))
        T[0] = [0, 0, 1]
        
        for i in range(1, 5):
            radius_vector = np.array([x[i]-x[0], y[i]-y[0], z[i]-z[0]])
            correcting_factor = self.curvature_correction(
                coplanarity_angle, radius_vector, n[:,0], n[:,i]
            )
            
            T[i, 0] = np.dot(radius_vector, l) * correcting_factor
            T[i, 1] = np.dot(radius_vector, p) * correcting_factor
            T[i, 2] = 1
        
        try:
            solution, _, _, _ = lstsq(T, gama)
            ql = solution[0]  # A coefficient
            qp = solution[1]  # B coefficient
            v_err = 0
        except:
            ql = 0
            qp = 0
            v_err = 1
        return ql,qp,v_err
    
def quadratic_3(self, x, y, z, l, p, n, coplanarity_angle, gama):
        """
        Calculate induced velocities using quadratic interpolation with 4 points
        """
        T = np.zeros((4, 4))
        
        for i in range(1, 4):
            radius_vector = np.array([x[i]-x[0], y[i]-y[0], z[i]-z[0]])
            correcting_factor = self.curvature_correction(
                coplanarity_angle, radius_vector, n[:,0], n[:,i]
            )
            
            T[i, 0] = np.dot(radius_vector, l)
            T[i, 1] = np.dot(radius_vector, p)
            T[i, 2] = T[i, 0]**2 + T[i, 1]**2
            T[i, 3] = 1
        
        try:
            solution = solve(T, gama)
            ql = solution[0]  # A coefficient
            qp = solution[1]  # B coefficient
            v_err = 0
        except:
            ql = 0
            qp = 0
            v_err = 1
        return ql,qp,v_err

def linear_3(self, points, l_vec, p_vec, normals, gama_nodal):
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
        T = np.zeros((4, 3))  # 转换矩阵
        T[0] = [0, 0, 1]      # 原点
        
        for i in range(1, 4):
            r_vec = np.array(points[i]) - np.array(points[0])
            coplanarity_angle = 0.1
            cf = self._curvature_correction(r_vec, normals[0], normals[i], coplanarity_angle)
            T[i,0] = np.dot(r_vec, l_vec) * cf  # l方向分量
            T[i,1] = np.dot(r_vec, p_vec) * cf  # p方向分量
            T[i,2] = 1.0                        # 常数项

        try:
            coeffs = solve(T, [gama_nodal[n-1] for n in points])
            return coeffs[0], coeffs[1], 0  # ql, qp, no error
        except:
            return 0, 0, 1  # 返回错误

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
        plane_unit_vector = np.cross(radius_unit_vector, average_normal)
        
        # Project normals to projection plane
        scalar = np.dot(normal_2, plane_unit_vector)
        normal_2_proj = normal_2 - scalar * plane_unit_vector
        scalar = np.dot(normal_1, plane_unit_vector)
        normal_1_proj = normal_1 - scalar * plane_unit_vector
        
        # Calculate angles of triangle
        ang1 = abs(np.arccos(np.dot(radius_unit_vector, normal_2_proj) / 
                          (np.linalg.norm(radius_unit_vector) * np.linalg.norm(normal_2_proj))) - pi/2)
        ang2 = abs(np.arccos(np.dot(radius_unit_vector, normal_1_proj) / 
                          (np.linalg.norm(radius_unit_vector) * np.linalg.norm(normal_1_proj))) - pi/2)
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