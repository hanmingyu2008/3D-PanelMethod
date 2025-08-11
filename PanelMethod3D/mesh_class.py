from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PanelMethod3D.plot_functions import set_axes_equal
import numpy as np
from PanelMethod3D.vector_class import Vector
from PanelMethod3D.panel_class import Panel, triPanel, quadPanel
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
        # self.shell_neighbours2 = self.locate_shells_adjacency2()
        
    
    @staticmethod
    def do_intersect(shell_i, shell_j): # 判断两个shell是否相交，用于求出邻接表(有公共边)
        
        return sum(node_id in shell_j for node_id in shell_i)>1    

    @staticmethod
    def do_intersect2(shell_i, shell_j): # 判断两个shell是否相交，用于求出邻接表(有公共顶点)
        
        return sum(node_id in shell_j for node_id in shell_i)>0  
    
    @staticmethod
    def do_involve(shell_i, node_id_list):
        
        return sum(node_id in shell_i for node_id in node_id_list) == len(node_id_list)
    
    def locate_shells_adjacency(self): 
        # 求一个网格的邻接表(有公共边)
        # 我使邻居列表里面的东西都是按照顺序来排列的
        shells = self.shells 
        neighbours = [] 
        for i, shell_i in enumerate(shells):
            neigh_old = []
            for j, shell_j in enumerate(shells): 
                if i != j and self.do_intersect(shell_i, shell_j): 
                    neigh_old.append(j) 
            neigh = []
            for k in range(len(shell_i)):
                k_next = (k+1)%len(shell_i)
                for j in neigh_old:
                    if self.do_involve(shells[j],[shell_i[k],shells[i][k_next]]):
                        neigh.append(j)
            neighbours.append(neigh)
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

    def refinement(self):
        nodes,shells = self.nodes,self.shells
        newnodes,shells = nodes.copy(),[]
        dic_edge = {}
        for she in shells:
            for k in range(len(she)):
                k_next = (k+1)%len(she)
                if frozenset({she[k],she[k_next]}) not in dic_edge:
                    dic_edge[frozenset({she[k],she[k_next]})] = len(newnodes)
                    x1,y1,z1 = nodes[she[k]]
                    x2,y2,z2 = nodes[she[k_next]]
                    newnodes.append(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2))
            if len(she) == 3:
                shells.append([she[0],dic_edge[frozenset({she[0],she[1]})],dic_edge[frozenset({she[0],she[2]})]])
                shells.append([she[1],dic_edge[frozenset({she[1],she[2]})],dic_edge[frozenset({she[1],she[0]})]])
                shells.append([she[2],dic_edge[frozenset({she[2],she[0]})],dic_edge[frozenset({she[2],she[1]})]])
            if len(she) == 4:
                n = len(newnodes)
                x1,y1,z1 = nodes[she[0]]
                x2,y2,z2 = nodes[she[1]]
                x3,y3,z3 = nodes[she[2]]
                x4,y4,z4 = nodes[she[3]]
                newnodes.append(((x1+x2+x3+x4)/4,(y1+y2+y3+y4)/4,(z1+z2+z3+z4)/4))
                shells.append([she[0],dic_edge[frozenset({she[0],she[1]})],n,dic_edge[frozenset({she[0],she[3]})]])
                shells.append([she[1],dic_edge[frozenset({she[1],she[2]})],n,dic_edge[frozenset({she[1],she[0]})]])
                shells.append([she[2],dic_edge[frozenset({she[2],she[3]})],n,dic_edge[frozenset({she[2],she[1]})]])
                shells.append([she[3],dic_edge[frozenset({she[3],she[0]})],n,dic_edge[frozenset({she[3],she[2]})]])
        return Mesh(newnodes,shells)
                          
class PanelMesh(Mesh):
    def __init__(self, nodes:list, shells:list):
        super().__init__(nodes, shells)
        self.panels = None
        self.panels_num = None
        self.panel_neighbours = self.shell_neighbours
        # self.panel_neighbours2 = self.shell_neighbours2
        # self.non_zero_panel = []
        self.CreatePanels()         
        
    def CreatePanels(self):
        panels = []
        for shell_id, shell in enumerate(self.shells):
            vertex = []
            for node_id in shell:
                node = self.nodes[node_id]
                vertex.append(Vector(node))

            if len(vertex) == 3:
                panel = triPanel(vertex[0],vertex[1],vertex[2])
                panels.append(panel)
                
            elif len(vertex) == 4:
                panel = quadPanel(vertex[0],vertex[1],vertex[2],vertex[3])
                panels.append(panel)
                
                
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
    
    def neighbour_of_neighbour(self,panel,i):
        # 我们要找到panel的第i号邻居的所有邻居当中,与panel相对的那一个
        nei = self.give_neighbours(panel)[i]
        if nei.num_vertices != 4:
            raise Exception("邻居必须是四边形才可以")
        k = self.give_neighbours(nei).index(panel)
        l = (k+2) % 4
        return self.give_neighbours(nei)[l]
    
    def refinement(self):
        # 加密,注意加密的时候也需要小心不要破坏了网格的规范(逆时针排列)
        nodes,shells = self.nodes,self.shells
        newnodes,shells = nodes.copy(),[]
        dic_edge = {}
        for she in self.shells:
            for k in range(len(she)):
                k_next = (k+1)%len(she)
                if frozenset({she[k],she[k_next]}) not in dic_edge:
                    dic_edge[frozenset({she[k],she[k_next]})] = len(newnodes)
                    x1,y1,z1 = nodes[she[k]]
                    x2,y2,z2 = nodes[she[k_next]]
                    newnodes.append(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2))
            if len(she) == 3:
                shells.append([she[0],dic_edge[frozenset({she[0],she[1]})],dic_edge[frozenset({she[0],she[2]})]])
                shells.append([she[1],dic_edge[frozenset({she[1],she[2]})],dic_edge[frozenset({she[1],she[0]})]])
                shells.append([she[2],dic_edge[frozenset({she[2],she[0]})],dic_edge[frozenset({she[2],she[1]})]])
                shells.append([dic_edge[frozenset({she[0],she[1]})],dic_edge[frozenset({she[1],she[2]})],
                               dic_edge[frozenset({she[2],she[0]})]])
            if len(she) == 4:
                n = len(newnodes)
                x1,y1,z1 = nodes[she[0]]
                x2,y2,z2 = nodes[she[1]]
                x3,y3,z3 = nodes[she[2]]
                x4,y4,z4 = nodes[she[3]]
                newnodes.append(((x1+x2+x3+x4)/4,(y1+y2+y3+y4)/4,(z1+z2+z3+z4)/4))
                shells.append([she[0],dic_edge[frozenset({she[0],she[1]})],n,dic_edge[frozenset({she[0],she[3]})]])
                shells.append([she[1],dic_edge[frozenset({she[1],she[2]})],n,dic_edge[frozenset({she[1],she[0]})]])
                shells.append([she[2],dic_edge[frozenset({she[2],she[3]})],n,dic_edge[frozenset({she[2],she[1]})]])
                shells.append([she[3],dic_edge[frozenset({she[3],she[0]})],n,dic_edge[frozenset({she[3],she[2]})]])
        return PanelMesh(newnodes,shells)
    
    def save(self,name):
        with open(name+".pan","w") as filee:
            print("# id starts from 0",file=filee)
            print("# nodes",file=filee)
            for x,y,z in self.nodes:
                print("n",x,y,z,file=filee)
            print("# panels",file=filee)
            for id,shell in enumerate(self.shells):
                if id in self.non_zero_panel:
                    if len(shell) == 3:
                        a,b,c = shell
                        print("f3",a,b,c,file=filee)
                    elif len(shell) == 4:
                        a,b,c,d = shell
                        print("f4",a,b,c,d,file=filee)
                    else:
                        raise Exception("len(shell) not in [3,4] !!!")
                else:
                    if len(shell) == 3:
                        a,b,c = shell
                        print("zf3",a,b,c,file=filee)
                    elif len(shell) == 4:
                        a,b,c,d = shell
                        print("zf4",a,b,c,d,file=filee)
                    else:
                        raise Exception("len(shell) not in [3,4] !!!")
            print("# neighbours",file=filee)
            for id,shell in enumerate(self.shells):
                lst = self.shell_neighbours[id]
                print("n"," ".join([str(n) for n in lst]))
    
class AeroMesh(Mesh):
    ## 现在我们希望可以自动判断trailing edge和sunction side等等
    ## 只需要输入哪些panel是wake,那些是body
    
    def __init__(self, nodes: list, shells: list, shells_ids:dict):
        super().__init__(nodes, shells)
        
        self.nodes = nodes
        self.shells = shells
        self.shells_ids = shells_ids
        '''
        self.TrailingEdge = {}
        self.wake_sheddingShells = {}
        
        self.set_shells_ids()
        self.find_TrailingEdge()
        self.set_wake_sheddingShells()
        '''
        self.find_TrailingEdge2()
        self.free_TrailingEdge()
        # self.free_LeadingEdge()
        # self.eliminate_main_surface_wing_tips_adjacency()
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
    
    def find_TrailingEdge2(self):
        SS_TE = []
        PS_TE = []
        self.wake_sheddingShells = {}
        temp = {}
        horseshore_tail = {}
        for wake_id in self.shells_ids["wake"]:
            temp[wake_id] = False
        num_temp = len(self.shells_ids["wake"])
        for wake_id in self.shells_ids["wake"]:
            lst = []
            for body_id in self.shells_ids["body"]:
                if self.do_intersect(self.shells[wake_id],self.shells[body_id]):
                    lst.append(body_id)
            if len(lst) == 1:
                print(wake_id)
                print(self.shells[wake_id])
                print(lst)
                print(self.shells[lst[0]])
                raise Exception("len(lst) == 1,为什么?!")
            if len(lst) == 0:
                continue
            if len(lst) > 2:
                raise Exception("len(lst) > 2,为什么?!")
            n1,n2 = lst
            self.wake_sheddingShells[n1] = [wake_id]
            horseshore_tail[n1] = self.shells[wake_id][-2:]
            self.wake_sheddingShells[n2] = [wake_id]
            horseshore_tail[n2] = self.shells[wake_id][-2:]
            temp[wake_id] = True
            num_temp -= 1
            v1 = Panel([Vector(self.nodes[i]) for i in self.shells[n1]]).r_cp 
            v2 = Panel([Vector(self.nodes[i]) for i in self.shells[n2]]).r_cp 
            wake_panel = Panel([Vector(self.nodes[t]) for t in self.shells[wake_id]])
            z1 = (v1-wake_panel.r_cp) * wake_panel.n
            z2 = (v2-wake_panel.r_cp) * wake_panel.n
            if z1 > z2:
                SS_TE.append(n1)
                PS_TE.append(n2)
            elif z1 < z2:
                SS_TE.append(n2)
                PS_TE.append(n1)
            else:
                print(wake_id)
                print("z1=z2,怎么会呢?")
        self.TrailingEdge = {"suction side": SS_TE, "pressure side": PS_TE}
        num_temp1 = num_temp
        num_temp2 = num_temp
        while num_temp1 != 0 or num_temp2 != 0:
            changed1 = False
            changed2 = False
            for wake_id in temp:
                if not temp[wake_id]:
                    a,b = self.shells[wake_id][:2]
                    for n in SS_TE:
                        if horseshore_tail[n] == [b,a]:
                            self.wake_sheddingShells[n].append(wake_id)
                            horseshore_tail[n] = self.shells[wake_id][-2:]
                            changed1 = True
                            num_temp1 -= 1
                    for n in PS_TE:
                        if horseshore_tail[n] == [b,a]:
                            self.wake_sheddingShells[n].append(wake_id)
                            horseshore_tail[n] = self.shells[wake_id][-2:]
                            changed2 = True
                            num_temp2 -= 1


            if not changed1 and not changed2:
                raise ValueError("wake不呈马鬃型")
            

    def free_TrailingEdge(self):
        """
        if trailing edge shells (or panels) share trailing edge nodes then
        the method "find_shell_neighbours" will assume that suction side trailing shells and 
        pressure side trailing shells are neighbours.
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
    
    
    def set_shells_ids(self):
        
        suction_side_shells_ids = self.find_suction_side()
        pressure_side_shells_ids = self.find_pressure_side()
        body_shells_ids = suction_side_shells_ids + pressure_side_shells_ids
        
        wake_shells_ids = self.shells_ids["wake"]
        
        self.shells_ids = {
            "body": body_shells_ids,
            "suction side": suction_side_shells_ids,
            "pressure side": pressure_side_shells_ids,
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
    
    def eliminate_body_wake_adjacency(self):
        # 做了点修改，以前"body"的位置是"main surface"
        if "body" in self.shells_ids and "wake" in self.shells_ids:
            self.eliminate_adjacency(
                self.shells_ids["body"], self.shells_ids["wake"]
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
    
    
    def locate_VSAERO_adjacency(self):
        
        """
        this function locates the neighbours of shells, following the notation of NASA Contractor Report 4023 
        "Program VSAERO theory Document,
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
        this function returns an array Ny X Nx array of panel ids,
         where Nx is the number of chordwise panels and Ny the number of spanwise panels.
        
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
    
    def __init__(self, nodes: list, shells: list, shells_ids: dict):
        super().__init__(nodes, shells, shells_ids)
        
        self.panels_ids = self.shells_ids
        self.wake_sheddingPanels = self.wake_sheddingShells
    
    def free_TrailingEdge(self):
        super().free_TrailingEdge()
        self.panel_neighbours = self.shell_neighbours

    def refinement(self):
        # 加密,注意加密的时候也需要小心不要破坏了网格的规范(逆时针排列)
        # 以及wake的规范更需要注意:       1 --- 0
        #                               |     |
        #                               2 --- 3
        nodes,shells,shells_ids = self.nodes,self.shells,self.shells_ids
        newnodes,shells,shells_ids_new = nodes.copy(),[],{"body":[],"wake":[]}
        dic_edge = {}
        for she_id in shells_ids["body"]:
            she = self.shells[she_id]
            for k in range(len(she)):
                k_next = (k+1)%len(she)
                if frozenset({she[k],she[k_next]}) not in dic_edge:
                    dic_edge[frozenset({she[k],she[k_next]})] = len(newnodes)
                    x1,y1,z1 = nodes[she[k]]
                    x2,y2,z2 = nodes[she[k_next]]
                    newnodes.append(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2))
        for she_id in shells_ids["wake"]:
            she = self.shells[she_id]
            assert(len(she) == 4)
            if frozenset({she[0],she[1]}) not in dic_edge:
                dic_edge[frozenset({she[0],she[1]})] = len(newnodes)
                x1,y1,z1 = nodes[she[0]]
                x2,y2,z2 = nodes[she[1]]
                newnodes.append(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2))
            if frozenset({she[2],she[3]}) not in dic_edge:
                dic_edge[frozenset({she[2],she[3]})] = len(newnodes)
                x1,y1,z1 = nodes[she[2]]
                x2,y2,z2 = nodes[she[3]]
                newnodes.append(((x1+x2)/2,(y1+y2)/2,(z1+z2)/2))

        for she_id in shells_ids["body"]:
            she = self.shells[she_id]
            if len(she) == 3:
                shells_ids_new["body"].append(len(shells))
                shells.append([she[0],dic_edge[frozenset({she[0],she[1]})],dic_edge[frozenset({she[0],she[2]})]])
                shells_ids_new["body"].append(len(shells))
                shells.append([she[1],dic_edge[frozenset({she[1],she[2]})],dic_edge[frozenset({she[1],she[0]})]])
                shells_ids_new["body"].append(len(shells))
                shells.append([she[2],dic_edge[frozenset({she[2],she[0]})],dic_edge[frozenset({she[2],she[1]})]])
                shells_ids_new["body"].append(len(shells))
                shells.append([dic_edge[frozenset({she[0],she[1]})],dic_edge[frozenset({she[1],she[2]})],
                               dic_edge[frozenset({she[2],she[0]})]])
            if len(she) == 4:
                n = len(newnodes)
                x1,y1,z1 = nodes[she[0]]
                x2,y2,z2 = nodes[she[1]]
                x3,y3,z3 = nodes[she[2]]
                x4,y4,z4 = nodes[she[3]]
                newnodes.append(((x1+x2+x3+x4)/4,(y1+y2+y3+y4)/4,(z1+z2+z3+z4)/4))
                shells_ids_new["body"].append(len(shells))
                shells.append([she[0],dic_edge[frozenset({she[0],she[1]})],n,dic_edge[frozenset({she[0],she[3]})]])
                shells_ids_new["body"].append(len(shells))
                shells.append([she[1],dic_edge[frozenset({she[1],she[2]})],n,dic_edge[frozenset({she[1],she[0]})]])
                shells_ids_new["body"].append(len(shells))
                shells.append([she[2],dic_edge[frozenset({she[2],she[3]})],n,dic_edge[frozenset({she[2],she[1]})]])
                shells_ids_new["body"].append(len(shells))
                shells.append([she[3],dic_edge[frozenset({she[3],she[0]})],n,dic_edge[frozenset({she[3],she[2]})]])

        for she_id in shells_ids["wake"]:
            she = self.shells[she_id]
            assert(len(she) == 4)
            shells_ids_new["wake"].append(len(shells))
            shells.append([dic_edge[frozenset({she[0],she[1]})],she[1],she[2],dic_edge[frozenset({she[2],she[3]})]])
            shells_ids_new["wake"].append(len(shells))
            shells.append([she[0],dic_edge[frozenset({she[0],she[1]})],dic_edge[frozenset({she[2],she[3]})],she[3]])

        return PanelAeroMesh(newnodes,shells,shells_ids_new)
    
    