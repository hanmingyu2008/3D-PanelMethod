from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plot_functions import set_axes_equal
import numpy as np
from vector_class import Vector
from panel_class import Panel, triPanel, quadPanel

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
       
if __name__=='__main__':
    from matplotlib import pyplot as plt
    from sphere import sphere
    nodes, shells = sphere(1, 10, 10, mesh_shell_type='triangular')
    sphere_mesh = PanelMesh(nodes, shells)    
    sphere_mesh.plot_panels()
    sphere_mesh.plot_mesh_bodyfixed_frame()
    
    