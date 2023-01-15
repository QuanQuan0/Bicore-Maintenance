# Core Maintenance
Graph formate: a .meta file contains the number of vertices in each sides and the number of edges; a .e file contains all the edges.

## Recompute
* Bi-core decomposition: ./abcore -ComShrDecom path_to_graph<br />path_to_graph: The path of the graph.<br />e.g. ./abcore -ComShrDecom ../data/example/

## Insertion maintenance
* To compute BiCore index: ./abcore -Build-BiCore path_to_graph<br />path_to_graph: The path of the graph.<br />e.g. ./abcore -Build-BiCore ../data/example/

* To compute bi-core numbers of each vertex: ./abcore -Build-BiCore-Number path_to_graph<br />path_to_graph: The path of the graph.<br />e.g. ./abcore -Build-BiCore-Number ../data/example/

* To insert edge with Edge-Insert: ./abcore -Edge-Insert path_to_graph vertex_1 vertex_2<br />path_to_graph: The path of the graph.<br />vertex_1/vertex_2: (vertex_1, vertex_2) is the inserted edge.<br />e.g. ./abcore -Edge-Insert ../data/example/ 1 1

## Deletion maintenance
* To compute BiCore index: ./abcore -Build-BiCore path_to_graph<br />path_to_graph: The path of the graph.<br />e.g. ./abcore -Build-BiCore ../data/example/

* To compute bi-core numbers of each vertex: ./abcore -Build-BiCore-Number path_to_graph<br />path_to_graph: The path of the graph.<br />e.g. ./abcore -Build-BiCore-Number ../data/example/

* To delete edge with Edge-Delete: ./abcore -Edge-Delete path_to_graph vertex_1 vertex_2<br />path_to_graph: The path of the graph.<br />vertex_1/vertex_2: (vertex_1, vertex_2) is the deleted edge.<br />e.g. ./abcore -Edge-Delete ../data/example/ 2 2