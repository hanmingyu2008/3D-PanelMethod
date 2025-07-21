from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import axes3d
from Algorithms import light_vector
from plot_functions import set_axes_equal, move_view
import numpy as np
from vector_class import Vector
from panel_class import Panel, triPanel, quadPanel

'''
网格Mesh类和PanelMesh、PanelAeroMesh类
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
    def do_intersect(shell_i, shell_j): # 判断两个shell是否相交，用于求出邻接表
        
        return sum(node_id in shell_j for node_id in shell_i)>1    

    @staticmethod
    def do_intersect2(shell_i, shell_j): # 判断两个shell是否相交，用于求出邻接表
        
        return sum(node_id in shell_j for node_id in shell_i)>0  
    
    def locate_shells_adjacency(self): # 求一个网格的邻接表 
        shells = self.shells 
        neighbours = [] 
        for i, shell_i in enumerate(shells):
            neighbours.append([]) 
            for j, shell_j in enumerate(shells): 
                if i != j and self.do_intersect(shell_i, shell_j): 
                    neighbours[-1].append(j) 
        
        return neighbours 
    
    def locate_shells_adjacency2(self): # 求一个网格的邻接表 
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
    
    def plot_shells(self, elevation=30, azimuth=-60):
        ax = plt.axes(projection='3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.view_init(elevation, azimuth)
        for shell in self.shells:
            x, y, z = [], [], []
            for id in shell:
                [x_i, y_i, z_i] = self.nodes[id]
                x.append(x_i)
                y.append(y_i)
                z.append(z_i)
            
            id = shell[0]
            [x_i, y_i, z_i] = self.nodes[id]
            x.append(x_i)
            y.append(y_i)
            z.append(z_i)    
            ax.plot3D(x, y, z, color='k')    
        
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()   
        
                          
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

    def give_neighbours(self, panel):
        
        neighbours_id_list = self.panel_neighbours[panel.id]
         
        neighbours_list = [self.panels[id] for id in neighbours_id_list]
        
        return neighbours_list
    
    def give_neighbours2(self, panel):
        
        neighbours_id_list = self.panel_neighbours2[panel.id]
         
        neighbours_list = [self.panels[id] for id in neighbours_id_list]
        
        return neighbours_list
    
    def give_neighbours3(self, panel, num):
        
        id_list = self.panel_neighbours[panel.id]

        id_list = sorted(id_list,key= lambda id: self.panels[id].n * panel.n,reverse=True)
         
        neighbours_list = [self.panels[id] for id in id_list[:num]]
        
        return neighbours_list
    
    def give_neighbours4(self, panel, num):
        
        id_list = self.panel_neighbours2[panel.id]

        id_list = sorted(id_list,key= lambda id: (self.panels[id].r_cp-panel.r_cp).norm(),reverse=False)
         
        neighbours_list = [self.panels[id] for id in id_list[:num]]
        
        return neighbours_list
    
    def plot_panels(self, elevation=30, azimuth=-60):
        ax = plt.axes(projection='3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.view_init(elevation, azimuth)
        
        for panel in self.panels:
            
            r_vertex = panel.r_vertex
            
            # plot panels
            if panel.num_vertices == 3:
                x = [r_vertex[0].x, r_vertex[1].x, r_vertex[2].x, r_vertex[0].x]
                y = [r_vertex[0].y, r_vertex[1].y, r_vertex[2].y, r_vertex[0].y]
                z = [r_vertex[0].z, r_vertex[1].z, r_vertex[2].z, r_vertex[0].z]
                ax.plot3D(x, y, z, color='k')
                
            elif panel.num_vertices == 4:
                
                x = [r_vertex[0].x, r_vertex[1].x, r_vertex[2].x, r_vertex[3].x,
                    r_vertex[0].x]
                y = [r_vertex[0].y, r_vertex[1].y, r_vertex[2].y, r_vertex[3].y,
                    r_vertex[0].y]
                z = [r_vertex[0].z, r_vertex[1].z, r_vertex[2].z, r_vertex[3].z,
                    r_vertex[0].z]
                ax.plot3D(x, y, z, color='k') 
                
            # plot normal vectors
            r_cp = panel.r_cp
            n = panel.n
            l = panel.l
            m = panel.m
            scale = 0.1 * panel.char_length
            n = n * scale
            l = l * scale
            m = m * scale
            ax.scatter(r_cp.x, r_cp.y, r_cp.z, color='k', s=5)
            ax.quiver(r_cp.x, r_cp.y, r_cp.z, n.x, n.y, n.z, color='r')
            ax.quiver(r_cp.x, r_cp.y, r_cp.z, l.x, l.y, l.z, color='b')
            ax.quiver(r_cp.x, r_cp.y, r_cp.z, m.x, m.y, m.z, color='g')
                
                   
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()

    def plot_mesh(self, elevation=30, azimuth=-60):
        shells = []
        vert_coords = []
        
        for panel in self.panels:
            shell=[]
            for r_vertex in panel.r_vertex:
                shell.append((r_vertex.x, r_vertex.y, r_vertex.z))
                vert_coords.append([r_vertex.x, r_vertex.y, r_vertex.z])
            shells.append(shell)
        
        
        light_vec = light_vector(magnitude=1, alpha=-45, beta=-45)
        face_normals = [panel.n for panel in self.panels]
        dot_prods = [-light_vec * face_normal for face_normal in face_normals]
        min = np.min(dot_prods)
        max = np.max(dot_prods)
        target_min = 0.2 # darker gray
        target_max = 0.6 # lighter gray
        shading = [(dot_prod - min)/(max - min) *(target_max - target_min) 
                   + target_min
                   for dot_prod in dot_prods]
        facecolor = plt.cm.gray(shading)
        
        ax = plt.axes(projection='3d')
        poly3 = Poly3DCollection(shells, facecolor=facecolor)
        ax.add_collection(poly3)
        ax.view_init(elevation, azimuth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        vert_coords = np.array(vert_coords)
        x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
        ax.set_xlim3d(x.min(), x.max())
        ax.set_ylim3d(y.min(), y.max())
        ax.set_zlim3d(z.min(), z.max())
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()

        
if __name__=='__main__':
    from matplotlib import pyplot as plt
    from sphere import sphere
    nodes, shells = sphere(1, 10, 10, mesh_shell_type='triangular')
    sphere_mesh = PanelMesh(nodes, shells)    
    sphere_mesh.plot_panels()
    sphere_mesh.plot_mesh_bodyfixed_frame()
    
    