from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plot_functions import set_axes_equal
import numpy as np
from vector_class import Vector
from panel_class import Panel, triPanel, quadPanel
from copy import deepcopy

'''
网格Mesh类和PanelMesh类
'''


class Mesh:
    
    def __init__(self, nodes:list, shells:list,
                 ):
        self.nodes = nodes
        self.shells = shells
        self.node_num = len(nodes)
        self.shell_num = len(shells)

        self.shell_neighbours = self.locate_shells_adjacency()
        self.shell_neighbours2 = self.locate_shells_adjacency2()
        
    
    @staticmethod
    def do_intersect(shell_i, shell_j): # 判断两个shell是否相交，用于求出邻接表(有公共边)
        
        return sum(node_id in shell_j for node_id in shell_i)>1    

    @staticmethod
    def do_intersect2(shell_i, shell_j): # 判断两个shell是否相交，用于求出邻接表(有公共顶点)
        
        return sum(node_id in shell_j for node_id in shell_i)>0  
    
    def locate_shells_adjacency(self): # 求一个网格的邻接表(有公共边)
        shells = self.shells 
        neighbours = [] 
        for i, shell_i in enumerate(shells):
            neighbours.append([]) 
            for j, shell_j in enumerate(shells): 
                if i != j and self.do_intersect(shell_i, shell_j): 
                    neighbours[-1].append(j) 
        
        return neighbours 
    
    def locate_shells_adjacency2(self): # 求一个网格的邻接表(有公共顶点)
        shells = self.shells 
        neighbours = [] 
        for i, shell_i in enumerate(shells):
            neighbours.append([]) 
            for j, shell_j in enumerate(shells): 
                if i != j and self.do_intersect2(shell_i, shell_j): 
                    neighbours[-1].append(j) 
        
        return neighbours

    def eliminate_adjacency(self, id_list1, id_list2): 
        
        for id in id_list1:
                        
            self.shell_neighbours[id] = [
                id for id in self.shell_neighbours[id] if id not in id_list2
            ]
                    
        for id in id_list2:
            
            self.shell_neighbours[id] = [
                id for id in self.shell_neighbours[id] if id not in id_list1
            ]
    
    def add_extra_neighbours(self):
         
        old_shell_neighbours = {}
        for id_i in range(self.shell_num):
            old_shell_neighbours[id_i] = self.shell_neighbours[id_i].copy()
            for id_j in self.shell_neighbours[id_i]:
                old_shell_neighbours[id_j] = self.shell_neighbours[id_j].copy()
        
        for id_i in range(self.shell_num):
            for id_j in old_shell_neighbours[id_i]:
                for id_k in old_shell_neighbours[id_j]: 
                    if id_k!=id_i and id_k not in self.shell_neighbours[id_i]:
                        self.shell_neighbours[id_i].append(id_k)   
                          
class PanelMesh(Mesh):
    def __init__(self, nodes:list, shells:list):
        super().__init__(nodes, shells)
        self.panels = None
        self.panels_num = None
        self.panel_neighbours = self.shell_neighbours
        self.panel_neighbours2 = self.shell_neighbours2
        self.CreatePanels()         
        
    def CreatePanels(self):
        panels = []
        for shell_id, shell in enumerate(self.shells):
            vertex = []
            for node_id in shell:
                node = self.nodes[node_id]
                vertex.append(Vector(node))

            if len(vertex) == 3:
                panels.append(triPanel(vertex[0], vertex[1], vertex[2]))
            elif len(vertex) == 4:
                panels.append(quadPanel(vertex[0], vertex[1],
                                        vertex[2], vertex[3]))
            
            panels[-1].id = shell_id
                
        self.panels = panels
        self.panels_num = len(panels)

    def give_neighbours(self, panel): # 返回某个panel的所有相邻panel(有相邻边)
        
        neighbours_id_list = self.panel_neighbours[panel.id]
         
        neighbours_list = [self.panels[id] for id in neighbours_id_list]
        
        return neighbours_list
    
    def give_neighbours2(self, panel): # 返回某个panel的所有相邻panel(有相邻顶点)
        
        neighbours_id_list = self.panel_neighbours2[panel.id]
         
        neighbours_list = [self.panels[id] for id in neighbours_id_list]
        
        return neighbours_list
    
    def give_neighbours3(self, panel, num): # 返回某个panel的所有相邻panel(有相邻边) 当中,法向量夹角更小的num个。用于后面求切向速度
        
        id_list = self.panel_neighbours[panel.id]

        id_list = sorted(id_list,key= lambda id: self.panels[id].n * panel.n,reverse=True)
         
        neighbours_list = [self.panels[id] for id in id_list[:num]]
        
        return neighbours_list
    
    def give_neighbours4(self, panel, num): # 返回某个panel的所有相邻panel(有相邻顶点) 当中,更近的num个。用于后面求切向速度。
        
        id_list = self.panel_neighbours2[panel.id]

        id_list = sorted(id_list,key= lambda id: (self.panels[id].r_cp-panel.r_cp).norm(),reverse=False)
         
        neighbours_list = [self.panels[id] for id in id_list[:num]]
        
        return neighbours_list
    
class AeroMesh(Mesh):
    
    def __init__(self, nodes: list, shells: list, nodes_ids:dict):
        super().__init__(nodes, shells)
        
        self.nodes_ids = nodes_ids
        self.shells_ids = {}
        self.TrailingEdge = {}
        self.wake_sheddingShells = {}
        
        self.set_shells_ids()
        self.find_TrailingEdge()
        self.set_wake_sheddingShells()
        
        self.free_TrailingEdge()
        # self.free_LeadingEdge()
        self.eliminate_main_surface_wing_tips_adjacency()
        self.eliminate_body_wake_adjacency()
    
    def find_TrailingEdge(self):
        te = self.nodes_ids["trailing edge"]
        ss = self.nodes_ids["suction side"]
        
        
        ps = self.nodes_ids["pressure side"]
        
        SS_TE = [shell_id for shell_id, shell in enumerate(self.shells)
                 if sum(node_id in te for node_id in shell)>1 
                 and sum(node_id in ss for node_id in shell) >2
                ]
        PS_TE = [shell_id for shell_id, shell in enumerate(self.shells)
                 if sum(node_id in te for node_id in shell)>1 
                 and sum(node_id in ps for node_id in shell) > 2
                ]
                
        self.TrailingEdge = {"suction side": SS_TE, "pressure side": PS_TE}

    def free_TrailingEdge(self):
        """
        if trailing edge shells (or panels) share trailing edge nodes then
        the method "find_shell_neighbours" will assume that suction side trailing shells and pressure side trailing shells are neighbours.
        In panel methods we assume that traling edge is a free edge.
        "free_TrailingEdge" method will remove false neighbour ids from the attribute shell_neighbour
        """
        
        list1 = [id for id in self.TrailingEdge["suction side"]]
        list2 = [id for id in self.TrailingEdge["pressure side"]]
        self.eliminate_adjacency(id_list1=list1, id_list2=list2)

    def free_LeadingEdge(self):
        num_shells_le = len(self.TrailingEdge["suction side"])
        LE_SS = self.shells_ids["suction side"][-num_shells_le:]
        LE_PS = self.shells_ids["pressure side"][0:num_shells_le]
        self.eliminate_adjacency(LE_SS, LE_PS)
        
    def find_suction_side(self):
        suction_side_nodes_ids = self.nodes_ids["suction side"]
        suction_side_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in suction_side_nodes_ids for node_id in shell)>2
        ]
        
        return suction_side_shells_ids
    
    def find_pressure_side(self):
        pressure_side_nodes_ids = self.nodes_ids["pressure side"]
        pressure_side_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in pressure_side_nodes_ids for node_id in shell)>2
        ]
        
        return pressure_side_shells_ids
    
    def find_right_wing_tip(self):
        right_wing_tip_nodes_ids = self.nodes_ids["right wing tip"]
        right_wing_tip_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in right_wing_tip_nodes_ids for node_id in shell)>2
        ]
        return right_wing_tip_shells_ids
    
    def find_left_wing_tip(self):
        left_wing_tip_nodes_ids =  self.nodes_ids["left wing tip"]
        left_wing_tip_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in left_wing_tip_nodes_ids for node_id in shell)>2
        ]
        return left_wing_tip_shells_ids
    
    def find_wake(self):
        wake_nodes_ids = self.nodes_ids["wake"]
        wake_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in wake_nodes_ids for node_id in shell)>2
        ]
        
        return wake_shells_ids
    
    def set_shells_ids(self):
        
        suction_side_shells_ids = self.find_suction_side()
        pressure_side_shells_ids = self.find_pressure_side()
        main_surface_shells_ids = suction_side_shells_ids\
            +pressure_side_shells_ids
        
        right_wing_tip_shells_ids = self.find_right_wing_tip()
        left_wing_tip_shells_ids = self.find_left_wing_tip()
        wing_tips_shells_ids = right_wing_tip_shells_ids \
            + left_wing_tip_shells_ids
        
        body_shells_ids = main_surface_shells_ids + wing_tips_shells_ids
        
        wake_shells_ids = self.find_wake()  
        
        self.shells_ids = {
            "body": body_shells_ids,
            "main surface": main_surface_shells_ids,
            "suction side": suction_side_shells_ids,
            "pressure side": pressure_side_shells_ids,
            "wing tips": wing_tips_shells_ids,
            "right tip": right_wing_tip_shells_ids,
            "left tip" : left_wing_tip_shells_ids,
            "wake": wake_shells_ids
        }
    
    def init_wake_sheddingShells(self):
        # if not empty dict
        if self.TrailingEdge:

            for j in range(len(self.TrailingEdge["suction side"])):

                self.wake_sheddingShells[self.TrailingEdge["suction side"][j]] = []

                self.wake_sheddingShells[self.TrailingEdge["pressure side"][j]] = []

    def set_wake_sheddingShells(self):
        
        self.init_wake_sheddingShells()
        
        for wake_line_id in range(len(self.nodes_ids["wake lines"])-1):
            
            line = self.nodes_ids["wake lines"][wake_line_id]
            next_line = self.nodes_ids["wake lines"][wake_line_id + 1]
            
            for te_shell_id in self.wake_sheddingShells:
                te_shell = self.shells[te_shell_id]
                
                if sum(
                    node_id in [line[0], next_line[0]]
                    for node_id in te_shell
                ) == 2:
                        
                        self.wake_sheddingShells[te_shell_id] = [
                            shell_id for shell_id in self.shells_ids["wake"]
                            if sum(
                                node_id in [*line, *next_line]
                                for node_id in self.shells[shell_id]
                            ) > 2
                        ]        
      
    def add_extra_neighbours(self):
         
        old_shell_neighbours = {}
        for id_i in self.shells_ids["body"]:
            old_shell_neighbours[id_i] = self.shell_neighbours[id_i].copy()
            for id_j in self.shell_neighbours[id_i]:
                old_shell_neighbours[id_j] = self.shell_neighbours[id_j].copy()
        
        for id_i in self.shells_ids["body"]:
            for id_j in old_shell_neighbours[id_i]:
                for id_k in old_shell_neighbours[id_j]: 
                    if id_k!=id_i and id_k not in self.shell_neighbours[id_i]:
                        self.shell_neighbours[id_i].append(id_k)      

    def eliminate_main_surface_wing_tips_adjacency(self):
        
        if "wing tips" in self.shells_ids and "main surface" in self.shells_ids:
            
            self.eliminate_adjacency(
                self.shells_ids["wing tips"], self.shells_ids["main surface"]
            )
    
    def eliminate_body_wake_adjacency(self):
        
        if "main surface" in self.shells_ids and "wake" in self.shells_ids:
            self.eliminate_adjacency(
                self.shells_ids["main surface"], self.shells_ids["wake"]
            )
    
    def give_near_root_shells_id(self):
        
        near_root_nodes_id = []
        for node_id, node in enumerate(self.nodes):
            if abs(node[1]) < 10**(-10):
                near_root_nodes_id.append(node_id)
        
        if self.shells_ids:
            body_shells_id = self.shells_ids["body"]
        else:
            body_shells_id = np.arange(len(self.shells))
        
        near_root_shells_id = []
        for shell_id in body_shells_id:
            for node_id in self.shells[shell_id]:
                if node_id in near_root_nodes_id:
                    near_root_shells_id.append(shell_id)
                    break
        
        return near_root_shells_id
    
    def give_leftSide_near_tip_shells_id(self):

        if self.shells_ids["left tip"]:
            left_wing_tip_shells_ids = self.shells_ids["left tip"]
        else:
            left_wing_tip_shells_ids = self.find_left_wing_tip()


        left_wing_tip_nodes_ids =  self.nodes_ids["left wing tip"]
        # wing tip shells' ids + near tip shells' ids
        left_side_near_tip_shells_id = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in left_wing_tip_nodes_ids for node_id in shell)>1
        ]

        # remove ids from wing tip shells
        left_side_near_tip_shells_id = [
            id for id in left_side_near_tip_shells_id if id not in left_wing_tip_shells_ids
        ]

        return left_side_near_tip_shells_id
    
    def locate_VSAERO_adjacency(self):
        
        """
        this function locates the neighbours of shells, following the notation of NASA Contractor Report 4023 "Program VSAERO theory Document,
        A Computer Program for Calculating Nonlinear Aerodynamic Characteristics
        of Arbitrary Configurations, Brian Maskew"
        """
        
        span_wise_nodes = len(self.nodes_ids["trailing edge"])
        chrod_wise_nodes = len(self.nodes_ids["main surface"])//span_wise_nodes
        span_wise_shells = span_wise_nodes - 1
        chrod_wise_shells = chrod_wise_nodes
        
        j_max = span_wise_shells - 1
        i_max = chrod_wise_shells - 1
        
        def shell_id(chord_wise_index, span_wise_index):            
            i = chord_wise_index
            j = span_wise_index
            shell_id = i*(j_max+1) + j
            return shell_id
              
              
        for i in range(1, i_max):
            for j in range(1, j_max):
                self.shell_neighbours[shell_id(i, j)] = [
                    shell_id(i, j-1),
                    shell_id(i+1, j),
                    shell_id(i, j+1),
                    shell_id(i-1, j)
                ]
        
        for i in range(1, i_max):
            j = 0
            self.shell_neighbours[shell_id(i, j)] = [
                -1,
                shell_id(i+1, j),
                shell_id(i, j+1),
                shell_id(i-1, j),
                
                shell_id(i, j+2)
            ]
            
            j = j_max
            self.shell_neighbours[shell_id(i, j)] = [
                shell_id(i, j-1),
                shell_id(i+1, j),
                -1,
                shell_id(i-1, j),
                
                shell_id(i, j-2)
            ]
            
        for j in range(1, j_max):
            i = 0
            self.shell_neighbours[shell_id(i, j)] = [
                shell_id(i, j-1),
                shell_id(i+1, j),
                shell_id(i, j+1),
                -1,
                
                shell_id(i+2, j)
            ]
            
            i = i_max
            self.shell_neighbours[shell_id(i, j)] = [
                shell_id(i, j-1),
                -1,
                shell_id(i, j+1),
                shell_id(i-1, j),
                
                shell_id(i-2, j)
            ]
            
        i, j = 0, 0
        self.shell_neighbours[shell_id(i, j)] = [
            -1,
            shell_id(i+1, j),
            shell_id(i, j+1),
            -1,
            
            shell_id(i, j+2),
            shell_id(i+2, j)
        ]
        
        i, j = 0, j_max
        self.shell_neighbours[shell_id(i, j)] = [
            shell_id(i, j-1),
            shell_id(i+1, j),
            -1,
            -1,
            
            shell_id(i, j-2),
            shell_id(i+2, j)
        ]
        
        i, j = i_max, 0
        self.shell_neighbours[shell_id(i, j)] = [
            -1,
            -1,
            shell_id(i, j+1),
            shell_id(i-1, j),
            
            shell_id(i, j+2),
            shell_id(i-2, j)
        ]
        
        i, j = i_max, j_max
        self.shell_neighbours[shell_id(i, j)] = [
            shell_id(i, j-1),
            -1,
            -1,
            shell_id(i-1, j),
            
            shell_id(i-2, j),
            shell_id(i, j-2)
        ]
      
    def give_ChordWiseStrips(self):
        """
        this function returns an array Ny X Nx array of panel ids, where Nx is the number of chordwise panels and Ny the number of spanwise panels.
        
        j-th row corresponds to j-th chordwise strip
        i-th column corresponds to i-th spanwise strip
        
        This function is meaningful only for quadrilateral structured meshes
        """

        span_wise_nodes = len(self.nodes_ids["trailing edge"])
        chrod_wise_nodes = len(self.nodes_ids["main surface"])//span_wise_nodes
        span_wise_shells = span_wise_nodes - 1
        chrod_wise_shells = chrod_wise_nodes

        j_max = span_wise_shells - 1
        i_max = chrod_wise_shells - 1

        def shell_id(chord_wise_index, span_wise_index):            
            i = chord_wise_index
            j = span_wise_index
            shell_id = i*(j_max+1) + j
            return shell_id

        ChordWiseStrips = np.zeros((j_max+1, i_max+1), dtype=int)
        for i in range(i_max + 1):
            for j in range(j_max + 1):
                ChordWiseStrips[j][i] = shell_id(i,j)

        return ChordWiseStrips
    
class PanelAeroMesh(AeroMesh, PanelMesh):
    
    def __init__(self, nodes: list, shells: list, nodes_ids: dict):
        super().__init__(nodes, shells, nodes_ids)
        
        self.panels_ids = self.shells_ids
        self.wake_sheddingPanels = self.wake_sheddingShells
    
    def free_TrailingEdge(self):
        super().free_TrailingEdge()
        self.panel_neighbours = self.shell_neighbours
      
    def give_near_root_panels(self):
        near_root_shells_id = super().give_near_root_shells_id()
        near_root_panels = []
        for id in near_root_shells_id:
            near_root_panels.append(self.panels[id])
        
        return near_root_panels
    
    def give_leftSide_near_root_panels(self):
        
        near_root_panels = self.give_near_root_panels()
        
        for panel in near_root_panels:
            if panel.r_cp.y < 0:
                near_root_panels.remove(panel)
        
        return near_root_panels

    def give_leftSide_near_tip_panels(self):

        leftSide_near_tip_panels = [
            self.panels[id] for id in self.give_leftSide_near_tip_shells_id()
        ]
        return leftSide_near_tip_panels
       
if __name__=='__main__':
    from matplotlib import pyplot as plt
    from sphere import sphere
    nodes, shells = sphere(1, 10, 10, mesh_shell_type='triangular')
    sphere_mesh = PanelMesh(nodes, shells)    
    sphere_mesh.plot_panels()
    sphere_mesh.plot_mesh_bodyfixed_frame()
    
    