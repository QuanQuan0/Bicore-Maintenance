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
    generated_edges = set()  
    attempts = 0  # 用于限制尝试次数

    while len(generated_edges) < num_edges and attempts < 10000:
        node = random.choice(list(B.nodes()))
        
        neighbors = list(B.neighbors(node))
        if neighbors:
            neighbor = random.choice(neighbors)
            edge = (node, neighbor)
            generated_edges.add(edge)

        attempts += 1

    return list(generated_edges)

def save_edges_to_file(edges, output_file_path):
    with open(output_file_path, 'w') as f:
        for edge in edges:
            f.write(f"{edge[0]} {edge[1]}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Input: python random_edges_generator.py <bipartite_graph_file> <num_edges> <output_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    num_edges_to_generate = int(sys.argv[2])
    output_file_path = sys.argv[3]

    B = load_bipartite_graph(file_path)
    random_edges = generate_random_edges(B, num_edges_to_generate)

    # 保存生成的边到输出文件
    save_edges_to_file(random_edges, output_file_path)

    print(f"生成的边已保存到 {output_file_path}")
