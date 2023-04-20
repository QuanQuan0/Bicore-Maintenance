#pragma once
#ifndef __BIGRAPH_H
#define __BIGRAPH_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <unordered_set>
#include "utility.h"

struct bicore_index_block {
	std::vector<vid_t> nodeset;
	bicore_index_block* next = NULL;
};

//struct bicore_index_block {
//	std::unordered_set<vid_t> nodeset;
//	bicore_index_block* next = NULL;
//};

struct skyline_block {
	std::vector<int> vector_set;
	//bicore_index_block* next = NULL;
};

struct bicore_index_block_dual_pointer {
	std::unordered_set<vid_t> nodeset;
	bicore_index_block_dual_pointer* horizontal_pointer = NULL;
	bicore_index_block_dual_pointer* vertical_pointer = NULL;
};

class Edge
{
public:
	Edge(int u_, int v_) { u = u_; v = v_; }
	bool operator<(const Edge &other) const
	{
		if (u == other.u)
			return v < other.v;
		return u < other.u;
	}

	int u;
	int v;
};


class DegreeNode
{
public:
	int id;
	int degree;
};

class BiGraph
{

public:

	BiGraph(std::string dir);
	BiGraph();
	~BiGraph() {}

	void addEdge(vid_t u, vid_t v);
	void deleteEdge(vid_t u, vid_t v);
	bool isEdge(vid_t u, vid_t v);
	num_t getV1Num() { return num_v1; }
	num_t getV2Num() { return num_v2; }
	num_t getV1Degree(vid_t u) { return degree_v1[u]; }
	num_t getV2Degree(vid_t u) { return degree_v2[u]; }
	std::vector<vid_t> & getV2Neighbors(vid_t u) { return neighbor_v2[u]; }
	std::vector<vid_t> & getV1Neighbors(vid_t u) { return neighbor_v1[u]; }
	void print();
	void print(bool hash);
	void printSum();
	void printCout();

public:

	void init(unsigned int num_v1, unsigned int num_v2);
	void loadGraph(std::string dir);
	void compressGraph();

	std::string dir;
	num_t num_v1;
	num_t num_v2;
	num_t num_edges;

	std::vector<std::vector<vid_t>> neighbor_v1;
	std::vector<std::vector<vid_t>> neighbor_v2;

	std::vector<std::unordered_set<vid_t>> neighborHash_v1;
	std::vector<std::unordered_set<vid_t>> neighborHash_v2;

	std::vector<int> degree_v1;
	std::vector<int> degree_v2;

	std::vector<num_t> core_v1;
	std::vector<num_t> core_v2;

	std::vector<std::map<std::vector<int>, int>> is_node_change_in_pair_v1;
	std::vector<std::map<std::vector<int>, int>> is_node_change_in_pair_v2;
	std::vector<int> store_skyline_index_u; std::vector<int> store_skyline_index_v;
    std::vector<int> store_skyline_index_u_beta; std::vector<int> store_skyline_index_v_beta;

public:

	//KKCore index left (x,*) right (*,x)
	std::vector<std::vector<int>> left_index;
	std::vector<std::vector<int>> right_index;
	std::vector<std::vector<bicore_index_block*>> g_bicore_index_u;
	std::vector<std::vector<bicore_index_block*>> g_bicore_index_v;
	int v1_max_degree;
	int v2_max_degree;
	std::vector<bool> left_delete;
	std::vector<bool> right_delete;
	std::vector<int> degree_of_left_node_initial;
    std::vector<int> degree_of_right_node_initial;
	// for dynamic update
	std::vector<std::vector<int>> left_index_old;
	std::vector<std::vector<int>> right_index_old;
	//BiGraph operator=(const BiGraph& g);
	int delta = -1;

public:
	int get_left_index_with_fixed_left_k(vid_t u, int left_k);
	//BiGraph& operator=(const BiGraph& g_);
};

extern void build_bicore_index(BiGraph&g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v);

extern void build_skyline_index(BiGraph&g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v);
//extern void build_skyline_index_check(BiGraph&g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u_check, std::vector<std::vector<bicore_index_block*>>& bicore_index_v_check, std::vector<skyline_block*>& skyline_index_u_check, std::vector<skyline_block*>& skyline_index_v_check,
//                         std::vector<skyline_block*>& skyline_index_v_reverse_check, std::vector<skyline_block*>& skyline_index_u_reverse_check);
extern void build_bicore_index_space_saver(BiGraph&g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v);

extern void build_bicore_index_space_saver_dual_pointer(BiGraph&g, std::vector<std::vector<bicore_index_block_dual_pointer*>>& bicore_index_u, std::vector<std::vector<bicore_index_block_dual_pointer*>>& bicore_index_v);

extern void retrieve_via_bicore_index(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	std::vector<bool>& left_node, std::vector<bool>& right_node, int alpha, int beta);

extern void retrieve_via_bicore_index_space_saver(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	std::vector<bool>& left_node, std::vector<bool>& right_node, int alpha, int beta);

extern void retrieve_via_bicore_index_space_saver_dual_pointer(BiGraph& g, std::vector<std::vector<bicore_index_block_dual_pointer*>>& bicore_index_u, std::vector<std::vector<bicore_index_block_dual_pointer*>>& bicore_index_v,
	std::vector<bool>& left_node, std::vector<bool>& right_node, int alpha, int beta);

extern void retrieve_via_bicore_index_inverse(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
	std::vector<bool>& left_node, std::vector<bool>& right_node, int alpha, int beta);

extern void inv(BiGraph& g);

extern void retrieve_bicore_via_bicore_number (BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v,
	std::vector<bool>& left_node, std::vector<bool>& right_node, int alpha, int beta);

extern void retrieve_node_via_index_skyline (BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v, batch_set& pair, vid_t node, bool is_left);

extern void retrieve_node_via_bicore_index (BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v
, batch_set& pair, vid_t node, bool is_left);

extern void retrieve_hierarchy_via_index_skyline (BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v, 
	std::vector<int>& left_node, std::vector<int>& right_node, std::queue<int> node_queue, std::queue<bool> is_left_queue, int alpha, int beta);

extern void retrieve_hierarchy_via_bicore_index(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v, 
	std::vector<bool>& left_node, std::vector<bool>& right_node, std::set<int>& left_set, std::set<int>& right_set, std::vector<int>& left_visit, std::vector<int>& right_visit, vid_t node, bool is_left, int alpha, int beta);

#endif  /* __BIGRAPH_H */
