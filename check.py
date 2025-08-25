import open3d as o3d
import numpy as np

def check_watertight(stl_file_path):
    """
    检查STL网格模型的水密性

    参数:
    stl_file_path (str): STL文件的路径

    返回:
    dict: 包含水密性、自交叠等检查结果的字典
    """

    # 读取STL文件
    mesh = o3d.io.read_triangle_mesh(stl_file_path)
    
    # 检查是否成功读取网格
    if mesh.is_empty():
        raise ValueError("无法读取网格或网格为空。")
    
    # 可视化网格（可选）
    # o3d.visualization.draw_geometries([mesh], window_name="Original Mesh")
    
    # 执行关键的拓扑检查

    print("="*50)
    print(f"网格水密性检查报告: {stl_file_path}")
    print("="*50)
    print(f"顶点数: {len(mesh.vertices)}")
    print(f"面片数: {len(mesh.triangles)}")

    is_watertight = mesh.is_watertight()

    print(f"是否为水密网格 (is_watertight): {'✅ 是' if is_watertight else '❌ 否'}")

    is_self_intersecting = mesh.is_self_intersecting()

    print(f"是否存在自交叠 (is_self_intersecting): {'❌ 是' if is_self_intersecting else '✅ 否'}")

    is_vertex_manifold = mesh.is_vertex_manifold()
    is_edge_manifold = mesh.is_edge_manifold(allow_boundary_edges=False)
    
    print(f"顶点流形检查 (is_vertex_manifold): {'✅ 是' if is_vertex_manifold else '❌ 否'}")
    print(f"边流形检查 (is_edge_manifold): {'✅ 是' if is_edge_manifold else '❌ 否'}")
    
    # 水密性综合判断
    if is_watertight and not is_self_intersecting:
        print("\n🎉 结论: 网格是水密的且无自交叠，可用于3D打印等工作。")
    else:
        print("\n⚠️  结论: 网格不是水密的或存在自交叠，需要修复。")
        # 如果不是水密的，检查有多少边界边
        if not is_watertight:
            # 注意：is_edge_manifold(allow_boundary_edges=True) 会返回边界边的数量
            _, boundary_edges = mesh.is_edge_manifold(allow_boundary_edges=True)
            print(f"   非流形边或边界边的数量: {boundary_edges}")

    # 将结果打包返回，方便程序后续处理
    result = {
        'is_watertight': is_watertight,
        'is_self_intersecting': is_self_intersecting,
        'is_vertex_manifold': is_vertex_manifold,
        'is_edge_manifold': is_edge_manifold,
        'num_vertices': len(mesh.vertices),
        'num_triangles': len(mesh.triangles),
    }
    return result

# 使用示例
if __name__ == "__main__":
    # 替换为你的STL文件路径
    file_path = "simplify_merged_surfaces (1).ply"
    
    try:
        check_watertight(file_path)
    except Exception as e:
        print(f"发生错误: {e}")