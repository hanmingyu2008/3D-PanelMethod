from collections import deque
from collections import defaultdict
import itertools

def is_homo_orient_neighbours(list1,list2):
    if sum(id in list1 for id in list2) < 2:
        raise Exception("list1与list2并不相邻啊?")
    elif sum(id in list1 for id in list2) > 2:
        return False
    for i in range(len(list1)):
        for j in range(len(list2)):
            if list1[i] == list2[j]:
                if list1[(i+1)%len(list1)] == list2[(j-1)%len(list2)]:
                    return True
                if list1[(i-1)%len(list1)] == list2[(j+1)%len(list2)]:
                    return True
            else:
                continue
    return False

def find_positive_faces(nodes, shells):
    """
    通过广度优先搜索（BFS）找到一个封闭三维几何表面中的所有正向面。

    Args:
        nodes: 一个包含三维顶点坐标的列表，形式为 [(x,y,z), ...]。
        shells: 一个包含面定义的列表，每个面由其顶点在nodes中的索引组成，
                形式为 [a,b,c] 或 [a,b,c,d]。

    Returns:
        一个包含所有正向面在shells列表中的索引的列表。
        如果找不到初始正向面，则会抛出ValueError。
    """

    if not shells:
        return []

    # Step 1: 标准化并统计每个面的出现次数
    # 标准化：对于一个面，总是从其最小索引的顶点开始，以确保循环等价的面被视为同一个。
    # 例如，[2,3,4,1] 和 [1,2,3,4] 都被标准化为 [1,2,3,4]
    normalized_shells_ids = []
    normalized_shells = []
    face_counts = defaultdict(int)
    for i,shell in enumerate(shells):
        # 找到最小的索引，并从那里开始重新排列顶点
        min_idx_val = min(shell)
        min_idx_pos = shell.index(min_idx_val)
        
        # 将列表旋转，使最小索引的顶点成为第一个
        # 这确保了 [1,2,3,4] 和 [2,3,4,1] 被视为同一个面
        normalized = shell[min_idx_pos:] + shell[:min_idx_pos]
        
        if normalized in normalized_shells:
            pass
        else:
            normalized_shells_ids.append(i)    
            normalized_shells.append(normalized)
            face_counts[tuple(sorted(normalized))] += 1
    
    # Step 2: 找到初始的正向面
    # 初始面被定义为只出现一次的面（既没有完全相同的副本，也没有反向副本）
    initial_face_index = -1
    for i in normalized_shells_ids:
        # 检查标准化后的面是否只出现一次
        if face_counts[tuple(sorted(shells[i]))] == 1:
            initial_face_index = i
            break
        
    if initial_face_index == -1:
        # raise ValueError("无法找到一个只出现一次的初始正向面。")
        pass
    # Step 3: 构建邻接列表
    # adj_list[i] 存储面 i 的所有邻居的索引
    adj_list = defaultdict(list)
    
    for i in normalized_shells_ids:
        shells_i = shells[i]
        for j in normalized_shells_ids:
            if i==j:
                continue
            shells_j = shells[j]
            if sum(id in shells_j for id in shells_i) > 1:
                adj_list[i].append(j)


    # Step 4: 执行广度优先搜索（BFS）
    initial_face_index = normalized_shells_ids[-1]
    queue = deque([initial_face_index])
    positive_faces_indices = {initial_face_index}  # 使用集合以获得 O(1) 的查找速度
    negative_faces_indices = set()

    while queue:
        current_face_idx = queue.popleft()
        current_face = shells[current_face_idx]

        for neighbor_idx in adj_list[current_face_idx]:
            if neighbor_idx not in positive_faces_indices and neighbor_idx not in negative_faces_indices:
                neighbor_face = shells[neighbor_idx]
                
                # 判断邻居是否为正向
                is_positive = is_homo_orient_neighbours(neighbor_face,current_face)
                
                if is_positive:
                    positive_faces_indices.add(neighbor_idx)
                    queue.append(neighbor_idx)
                else:
                    pass# negative_faces_indices.add(neighbor_idx)

    return sorted(list(positive_faces_indices))

def find_unique_faces(shells):
    new_shells = []
    ref = []
    for shell in shells:
        if sorted(shell) not in ref:
            new_shells.append(shell)
            ref.append(sorted(shell))
    return new_shells

def preprocessing_off(filename,outputname):
    with open(filename+".off","r") as filee:
        lines = []
        for line in filee:
            lines.append(line)
        a,b,c = [int(x) for x in lines[1].split()]
        shells = []
        for id in range(a+2,a+b+2):
            shells.append([int(x) for x in lines[id].split()[1:]])
        new_shells = find_unique_faces(shells)
    with open(outputname+".off","w") as filee:
        print("OFF",file=filee)
        print(a,len(new_shells),c,file=filee)
        for id in range(2,a+2):
            print(lines[id],end="",file=filee)
        for shell in new_shells:
            shell.insert(0,len(shell))
            print(" ".join([str(x) for x in shell]),file=filee)

if __name__ == "__main__":

    # 示例用法
    nodes_example = [(0,0,0), (1,0,0), (1,1,0), (0,1,0),
                    (0,0,1), (1,0,1), (1,1,1), (0,1,1)]

    # 这是一个立方体的面，其中一个面是唯一的，
    # 另一些面是重复或反向的，以模拟你的场景。
    shells_example = [
        # 底部面 (正向)
        [0, 1, 2, 3],
        # 顶部面 (反向) - 应该是 [4,5,6,7]
        [7, 6, 5, 4],
        # 侧面 1
        [0, 4, 5, 1],
        # 侧面 1 的反向副本
        [1, 5, 4, 0],
        # 侧面 2
        [1, 5, 6, 2],
        # 侧面 3
        [2, 6, 7, 3],
        # 侧面 4
        [3, 7, 4, 0],
        # 顶部面的正向版本，但多了一个
        [4, 5, 6, 7],
        # 另外一个底部面，方向相反
        [3, 2, 1, 0],
        [3, 0, 1, 2],
        [6, 7, 8, 9],
        [6, 9, 8, 7]
    ]

    try:
        positive_faces = find_positive_faces(nodes_example, shells_example)
        print("找到的正向面的索引为:", positive_faces)
    except ValueError as e:
        print(e)