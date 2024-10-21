import networkx as nx
import random
import sys

# 读取二部图文件
def load_bipartite_graph(file_path):
    B = nx.Graph()
    with open(file_path, 'r') as f:
        for line in f:
            u, v = line.strip().split()
            B.add_edge(u, v)
    return B

# 随机生成边
def generate_random_edges(B, num_edges):
    generated_edges = set()  # 用于存储生成的边，避免重复
    attempts = 0  # 用于限制尝试次数

    while len(generated_edges) < num_edges and attempts < 10000:
        # 随机选择一个点
        node = random.choice(list(B.nodes()))
        
        # 获取该点的所有邻居
        neighbors = list(B.neighbors(node))

        # 如果点有邻居，随机选择一个邻居并生成边
        if neighbors:
            neighbor = random.choice(neighbors)
            edge = (node, neighbor)
            generated_edges.add(edge)

        attempts += 1

    return list(generated_edges)

# 将生成的边写入文件
def save_edges_to_file(edges, output_file_path):
    with open(output_file_path, 'w') as f:
        for edge in edges:
            f.write(f"{edge[0]} {edge[1]}\n")

# 从命令行参数获取文件路径和边的数量
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Input: python random_edges_generator.py <bipartite_graph_file> <num_edges> <output_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    num_edges_to_generate = int(sys.argv[2])
    output_file_path = sys.argv[3]

    B = load_bipartite_graph(file_path)

    # 生成指定数量的随机边
    random_edges = generate_random_edges(B, num_edges_to_generate)

    # 保存生成的边到输出文件
    save_edges_to_file(random_edges, output_file_path)

    print(f"生成的边已保存到 {output_file_path}")
