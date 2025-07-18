import numpy as np
from scipy.linalg import solve, lstsq

class VelocityCalculator:
    """
    速度计算器类 - 计算由双偶极子强度诱导的速度场
    核心方法：
        1. 节点插值法（Nodal Interpolation）
        2. 表面插值法（Surface Interpolation）
    """
    def __init__(self):
        self.velo_err = 0       # 错误标志：0=正常, 1=计算错误
        self.v = None           # 总速度矩阵 [面板数 × 工况数]
        self.vx = None          # x方向速度分量矩阵
        self.vy = None          # y方向速度分量矩阵
        self.vz = None          # z方向速度分量矩阵

    def calculate_velocities(self, params):
        """
        主计算函数 - 根据输入参数计算诱导速度场
        参数说明：
            params: 包含以下键的字典：
                'gama'              - 双偶极子强度矩阵 [面板数 × 工况数]
                'node_num'          - 总节点数
                'panel_num'         - 总面板数 
                'panel_num_no_wake' - 无尾流面板数
                'velorder'          - 插值阶数 (0=节点插值, 1=线性, 2=二次)
                'method'            - 计算方法 (0=源/双偶极子, 1=仅双偶极子)
                'coplanarity_angle' - 共面性角度阈值（弧度）
                ...（其他几何参数）
        """
        try:
            # === 初始化速度场 ===
            self.v = np.zeros((params['panel_num_no_wake'], params['case_num']))
            self.vx = np.zeros_like(self.v)
            self.vy = np.zeros_like(self.v)
            self.vz = np.zeros_like(self.v)

            self._surface_interpolation(params)

        except Exception as e:
            print(f"[ERROR] 速度计算失败: {str(e)}")
            self.velo_err = 1

    def _surface_interpolation(self, params):
        """表面插值法核心逻辑"""
        for i in range(params['panel_num_no_wake']):
            if params['panel_type'][i] in (20, 21):  # 跳过虚拟面板
                continue

            # 获取相邻面板数据
            neighbors = self._get_neighbor_panels(i, params)
            
            # 根据邻居数量选择插值方法
            if len(neighbors) == 4:
                if params['velorder'] == 2:
                    ql, qp, v_err = self._quadratic_5(neighbors, params)
                else:
                    ql, qp, v_err = self._linear_5(neighbors, params)
            elif len(neighbors) == 3:
                if params['velorder'] == 2:
                    ql, qp, v_err = self._quadratic_4(neighbors, params)
                else:
                    ql, qp, v_err = self._linear_4(neighbors, params)
            else:
                ql, qp, v_err = self._linear_3(neighbors, params)

            # 处理计算错误
            if v_err:
                ql, qp = self._handle_error(params['method'], params['speed'])

            # 转换到全局坐标系
            self._compute_global_velocity(i, params, ql, qp)

    def _compute_global_velocity(self, i, params, ql, qp):
        """将局部速度(ql,qp)转换为全局坐标系(vx,vy,vz)"""
        for j in range(params['case_num']):
            if params['method'] == 0:  # 源/双偶极子组合
                # 计算自由流分量
                gl = np.dot(params['l'][i], [params['speed_x'][j], params['speed_y'][j], params['speed_z'][j]])
                gp = np.dot(params['p'][i], [params['speed_x'][j], params['speed_y'][j], params['speed_z'][j]])
                
                # 计算扰动速度
                self.vx[i,j] = (gl - ql[j]) * params['l'][i][0] + (gp - qp[j]) * params['p'][i][0]
                self.vy[i,j] = (gl - ql[j]) * params['l'][i][1] + (gp - qp[j]) * params['p'][i][1]
                self.vz[i,j] = (gl - ql[j]) * params['l'][i][2] + (gp - qp[j]) * params['p'][i][2]
                
                # 总速度
                self.v[i,j] = np.sqrt(self.vx[i,j]**2 + self.vy[i,j]**2 + self.vz[i,j]**2)
                
            else:  # 仅双偶极子方法
                self.vx[i,j] = -ql[j] * params['l'][i][0] - qp[j] * params['p'][i][0]
                self.vy[i,j] = -ql[j] * params['l'][i][1] - qp[j] * params['p'][i][1]
                self.vz[i,j] = -ql[j] * params['l'][i][2] - qp[j] * params['p'][i][2]
                self.v[i,j] = np.sqrt(ql[j]**2 + qp[j]**2)

    # ====================== 插值核心函数 ======================

    def _quadratic_5(self, x, y, z, l, p, n, coplanarity_angle, gama):
        """
        Calculate induced velocities using quadratic interpolation with 5 points
        """
        T = np.zeros((5, 5))
        T[0] = [0, 0, 0, 0, 1]
        
        for i in range(1, 5):
            radius_vector = np.array([x[i]-x[0], y[i]-y[0], z[i]-z[0]])
            correcting_factor = self._curvature_correction(
                coplanarity_angle, radius_vector, n[:,0], n[:,i]
            )
            
            T[i, 0] = np.dot(radius_vector, l) * correcting_factor
            T[i, 1] = np.dot(radius_vector, p) * correcting_factor
            T[i, 2] = T[i, 0]**2
            T[i, 3] = T[i, 1]**2
            T[i, 4] = 1
        
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
    
    def _linear_5(self, x, y, z, l, p, n, coplanarity_angle, gama):
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
    
    def _quadratic_4(self, x, y, z, l, p, n, coplanarity_angle, gama):
        """
        Calculate induced velocities using quadratic interpolation with 4 points
        """
        T = np.zeros((4, 4))
        T[0] = [0, 0, 0, 1]
        
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

    def _linear_4(self, points, l_vec, p_vec, normals, gama_nodal):
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

    def _curvature_correction(self, coplanarity_angle, radius_vector, normal_1, normal_2):
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

    def _calculate_weights(self, weights, params):
        """计算逆距离权重 (需要根据实际需求实现)"""
        pass

    def _get_panel_coords(self, panel_idx, params):
        """获取面板节点坐标"""
        nodes = [
            params['node1'][panel_idx],
            params['node2'][panel_idx], 
            params['node3'][panel_idx]
        ]
        if params['panel_type'][panel_idx] in (1, 20):
            nodes.append(params['node4'][panel_idx])
        return [(params['x'][n-1], params['y'][n-1], params['z'][n-1]) for n in nodes]

    def _handle_error(self, method, speed):
        """处理插值错误"""
        if method == 0:
            return np.zeros_like(speed), np.zeros_like(speed)
        else:
            return np.sqrt(speed), np.sqrt(speed)