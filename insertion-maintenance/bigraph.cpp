#include "bigraph.h"
#include <string>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <list>
#include <vector>
#include <unordered_set>
#include <map>
#include <set>
#include <cstdlib>
#include "utility.h"
using namespace std;

#define LINE_LENGTH 1000



BiGraph::BiGraph(string dir)
{
	num_v1 = 0;
	num_v2 = 0;
	num_edges = 0;

	neighbor_v1.clear();
	neighbor_v2.clear();


	degree_v1.clear();
	degree_v2.clear();

	core_v1.clear();
	core_v2.clear();

	//KKCore index left (x,*) right (*,x)
	left_index.clear();
	right_index.clear();
	v1_max_degree = 0;
	v2_max_degree = 0;
	delta = -1;
	this->dir = dir;
	loadGraph(dir);
}

BiGraph::BiGraph() {
	dir = "";
	num_v1 = 0;
	num_v2 = 0;
	num_edges = 0;

	neighbor_v1.clear();
	neighbor_v2.clear();


	degree_v1.clear();
	degree_v2.clear();

	core_v1.clear();
	core_v2.clear();

	//KKCore index left (x,*) right (*,x)
	left_index.clear();
	right_index.clear();
	v1_max_degree = 0;
	v2_max_degree = 0;
	delta = -1;
}

void BiGraph::print()
{
	string bigraphE = dir + "graph.e";
	string bigraphMeta = dir + "graph.meta";

	FILE *graphEF = fopen(bigraphE.c_str(), "w");
	FILE *graphMetaF = fopen(bigraphMeta.c_str(), "w");

	fprintf(graphMetaF, "%d\n%d\n%d\n", num_v1, num_v2, num_edges);
	fclose(graphMetaF);
	for (int i = 0; i < num_v1; ++i)
	{
		for (int j = 0; j < neighbor_v1[i].size(); ++j)
		{
			fprintf(graphEF, "%d %d\n", i, neighbor_v1[i][j]);
		}
	}
	fclose(graphEF);
}

void BiGraph::print(bool hash) {
	string bigraphE = dir + "graph.e";
	string bigraphMeta = dir + "graph.meta";

	FILE *graphEF = fopen(bigraphE.c_str(), "w");
	FILE *graphMetaF = fopen(bigraphMeta.c_str(), "w");

	fprintf(graphMetaF, "%d\n%d\n%d\n", neighborHash_v1.size(), neighborHash_v2.size(), num_edges);
	fclose(graphMetaF);
	for (int i = 0; i < neighborHash_v1.size(); ++i)
	{
		for (auto j = neighborHash_v1[i].begin(); j != neighborHash_v1[i].end(); ++j)
		{
			fprintf(graphEF, "%d %d\n", i, *j);
		}
	}
	fclose(graphEF);
}

void BiGraph::printCout()
{
	cout << "\nBiGraph: " << endl;
	for (int i = 0; i < num_v1; ++i)
	{
		cout << i << ": ";
		if (neighbor_v1[i].size() == 0)
		{
			cout << "compress error" << endl;
			exit(1);
		}
		for (int j = 0; j < neighbor_v1[i].size(); ++j)
		{
			cout << neighbor_v1[i][j] << ", ";
		}
		cout << endl;
	}
	cout << endl;

}

void BiGraph::printSum()
{
	cout << "\nBiGraph Sum: " << endl;
	cout << "num_v1: " << num_v1 << endl;
	cout << "num_v2: " << num_v2 << endl;
	cout << "edge: " << num_edges << endl;

}

void BiGraph::init(unsigned int num1, unsigned int num2)
{
	num_v1 = num1;
	num_v2 = num2;
	num_edges = 0;

	neighbor_v1.resize(num_v1);
	neighbor_v2.resize(num_v2);

	degree_of_left_node_initial.resize(num_v1);
	degree_of_right_node_initial.resize(num_v2);

	is_node_change_in_pair_v1.resize(num_v1);
    is_node_change_in_pair_v2.resize(num_v2);
    store_skyline_index_u.resize(num_v1); store_skyline_index_v.resize(num_v2);
    store_skyline_index_u_beta.resize(num_v1); store_skyline_index_v_beta.resize(num_v2);

	fill_n(degree_of_left_node_initial.begin(), num_v1, 0);
	fill_n(degree_of_right_node_initial.begin(), num_v2, 0);

	degree_v1.resize(num_v1);
	degree_v2.resize(num_v2);

	fill_n(degree_v1.begin(), num_v1, 0);
	fill_n(degree_v2.begin(), num_v2, 0);

	left_delete.resize(num_v1);
	right_delete.resize(num_v2);
}

void BiGraph::loadGraph(string dir)
{
	unsigned int n1, n2;
	unsigned int edges = 0;
	int u, v;
	int r;

	string metaFile = dir + "graph.meta";
	string edgeFile = dir + "graph.e";

	FILE * metaGraph = fopen(metaFile.c_str(), "r");
	FILE * edgeGraph = fopen(edgeFile.c_str(), "r");

	if (fscanf(metaGraph, "%d\n%d", &n1, &n2) != 2)
	{
		fprintf(stderr, "Bad file format: n1 n2 incorrect\n");
		exit(1);
	}

	fprintf(stdout, "n1: %d, n2: %d\n", n1, n2);

	init(n1, n2);

	while ((r = fscanf(edgeGraph, "%d %d", &u, &v)) != EOF)
	{
		//fprintf(stderr, "%d, %d\n", u, v);
		if (r != 2)
		{
			fprintf(stderr, "Bad file format: u v incorrect\n");
			exit(1);
		}

		addEdge(u, v);
		//num_edges++;
	}

	fclose(metaGraph);
	fclose(edgeGraph);

	for (int i = 0; i < num_v1; ++i)
	{
		neighbor_v1[i].shrink_to_fit();
		sort(neighbor_v1[i].begin(), neighbor_v1[i].end());

	}
	for (int i = 0; i < num_v2; ++i)
	{
		neighbor_v2[i].shrink_to_fit();
		sort(neighbor_v2[i].begin(), neighbor_v2[i].end());
	}

	//neighborHash_v1.resize(num_v1);
	//neighborHash_v2.resize(num_v2);

	//	for (int i = 0; i < num_v1; ++i)
	//	{
	//		neighborHash_v1[i].reserve(10*neighbor_v1[i].size());
	//		for (int j = 0; j < neighbor_v1[i].size(); ++j)
	//		{
	//			int v = neighbor_v1[i][j];
	//			neighborHash_v1[i].insert(v);
	//		}
	//	}

	//	for (int i = 0; i < num_v2; ++i)
	//	{
	//		neighborHash_v2[i].reserve(2*neighbor_v2[i].size());
	//		for (int j = 0; j < neighbor_v2[i].size(); ++j)
	//		{
	//			int v = neighbor_v2[i][j];
	//			neighborHash_v2[i].insert(v);
	//		}
	//	}

}

void BiGraph::addEdge(vid_t u, vid_t v)
{

	neighbor_v1[u].push_back(v);
	++degree_v1[u];
	if (degree_v1[u] > v1_max_degree) v1_max_degree = degree_v1[u];
	neighbor_v2[v].push_back(u);
	++degree_v2[v];
	if (degree_v2[v] > v2_max_degree) v2_max_degree = degree_v2[v];
	num_edges++;
}

// not change max_degree
void BiGraph::deleteEdge(vid_t u, vid_t v)
{
	for (int i = 0; i < degree_v1[u]; ++i)
	{
		int vv = neighbor_v1[u][i];
		if (vv == v)
		{
			swap(neighbor_v1[u][i], neighbor_v1[u][degree_v1[u] - 1]);
			--degree_v1[u];
			neighbor_v1[u].pop_back();
			num_edges--;//only once!!!
			break;
		}
	}

	if (degree_v1[u] + 1 == v1_max_degree) {
		v1_max_degree = 0;
		for (auto d : degree_v1) {
			v1_max_degree = v1_max_degree < d ? d : v1_max_degree;
		}
	}

	for (int i = 0; i < degree_v2[v]; ++i)
	{
		int uu = neighbor_v2[v][i];
		if (uu == u)
		{
			swap(neighbor_v2[v][i], neighbor_v2[v][degree_v2[v] - 1]);
			--degree_v2[v];
			neighbor_v2[v].pop_back();
			break;
		}
	}

	if (degree_v2[v] + 1 == v2_max_degree) {
		v2_max_degree = 0;
		for (auto d : degree_v2) {
			v2_max_degree = v2_max_degree < d ? d : v2_max_degree;
		}
	}
}

bool BiGraph::isEdge(vid_t u, vid_t v)
{
	/*if (binary_search(neighbor_v1[u].begin(),
		neighbor_v1[u].begin() + degree_v1[u], v))
		return true;
	else
		return false;*/
	/*if (neighborHash_v1[u].find(v) == neighborHash_v1[u].end())
	{
	return false;
	}
	else
	return true;*/
	for (auto it = neighbor_v1[u].begin(); it != neighbor_v1[u].end(); it++) {
		if (*it == v) return true;
	}
	return false;
}

void BiGraph::compressGraph()
{
	vector<unordered_set<vid_t>> n_neighborHash_v1;
	vector<unordered_set<vid_t>> n_neighborHash_v2;

	n_neighborHash_v1.resize(num_v1);
	n_neighborHash_v2.resize(num_v2);

		for (int i = 0; i < num_v1; ++i)
		{
			n_neighborHash_v1[i].reserve(10*neighbor_v1[i].size());
			for (int j = 0; j < neighbor_v1[i].size(); ++j)
			{
				int v = neighbor_v1[i][j];
				n_neighborHash_v1[i].insert(v);
			}
		}

		for (int i = 0; i < num_v2; ++i)
		{
			n_neighborHash_v2[i].reserve(2*neighbor_v2[i].size());
			for (int j = 0; j < neighbor_v2[i].size(); ++j)
			{
				int v = neighbor_v2[i][j];
				n_neighborHash_v2[i].insert(v);
			}
		}

	swap(n_neighborHash_v1, neighborHash_v1);
	swap(n_neighborHash_v2, neighborHash_v2);

	for (int i = 0; i < num_v1; ++i)
	{
	if (neighbor_v1[i].size() != degree_v1[i])
	cout << "degree error" << endl;
	}

	for (int i = 0; i < num_v2; ++i)
	{
	if (neighbor_v2[i].size() != degree_v2[i])
	cout << "degree error" << endl;
	}

	cout<<"degree correct"<<endl;

}

int BiGraph::get_left_index_with_fixed_left_k(vid_t u, int left_k) {
	if (left_index[u].size() > left_k) return left_index[u][left_k];
	else return 0;
}

/*BiGraph& BiGraph::operator=(const BiGraph& g_) {
	num_v1 = g_.num_v1;
	num_v2 = g_.num_v2;
	num_edges = g_.num_edges;
	neighbor_v1.clear(); neighbor_v1.resize(num_v1);
	neighbor_v2.clear(); neighbor_v2.resize(num_v2);
	degree_v1.clear(); degree_v1.resize(num_v1);
	degree_v2.clear(); degree_v2.resize(num_v2);
	left_delete.clear(); left_delete.resize(num_v1);
	right_delete.clear(); right_delete.resize(num_v2);
	for (vid_t u = 0; u < num_v1; u++) {
		degree_v1[u] = g_.degree_v1[u];
		left_delete[u] = g_.left_delete[u];
		neighbor_v1[u].resize(g_.neighbor_v1[u].size());
		for (int i = 0; i < neighbor_v1[u].size(); i++) {
			neighbor_v1[u][i] = g_.neighbor_v1[u][i];
		}
	}
	for (vid_t v = 0; v < num_v2; v++) {
		degree_v2[v] = g_.degree_v2[v];
		right_delete[v] = g_.right_delete[v];
		neighbor_v2[v].resize(g_.neighbor_v2[v].size());
		for (int i = 0; i < neighbor_v2[v].size(); i++) {
			neighbor_v2[v][i] = g_.neighbor_v2[v][i];
		}
	}
	v1_max_degree = g_.v1_max_degree;
	v2_max_degree = g_.v2_max_degree;
	return *this;
}*/

void build_bicore_index(BiGraph&g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v) {
	bicore_index_u.clear(); bicore_index_u.resize(g.v1_max_degree + 1); bicore_index_v.clear(); bicore_index_v.resize(g.v2_max_degree + 1);
	// build left
	vector<int> beta_m; beta_m.resize(g.v1_max_degree + 1);
	for (vid_t u = 0; u < g.num_v1; u++) {
		for (int alpha = 1; alpha < g.left_index[u].size(); alpha++) {
			int beta = g.left_index[u][alpha];
			if (beta_m[alpha] < beta) beta_m[alpha] = beta;
		}
	}
	for (int alpha = 1; alpha < beta_m.size(); alpha++) {
		bicore_index_u[alpha].resize(beta_m[alpha] + 1);
		for (int kkk = 1; kkk < bicore_index_u[alpha].size(); kkk++) {
			bicore_index_u[alpha][kkk] = new bicore_index_block;
		}
	}
	for (vid_t u = 0; u < g.num_v1; u++) {
		for (int alpha = 1; alpha < g.left_index[u].size(); alpha++) {
			int beta = g.left_index[u][alpha];
			bicore_index_u[alpha][beta]->nodeset.push_back(u);
		}
	}
	for (int alpha = 1; alpha < bicore_index_u.size(); alpha++) {
		bicore_index_block* pre = NULL;
		for (int beta = bicore_index_u[alpha].size() - 1; beta > 0; beta--) {
			bicore_index_u[alpha][beta]->next = pre;
			if (bicore_index_u[alpha][beta]->nodeset.size() > 0) pre = bicore_index_u[alpha][beta];
		}
	}
	// build right
	vector<int> alpha_m; alpha_m.resize(g.v2_max_degree + 1);
	for (vid_t v = 0; v < g.num_v2; v++) {
		for (int beta = 1; beta < g.right_index[v].size(); beta++) {
			int alpha = g.right_index[v][beta];
			if (alpha_m[beta] < alpha) alpha_m[beta] = alpha;
		}
	}
	for (int beta = 1; beta < alpha_m.size(); beta++) {
		bicore_index_v[beta].resize(alpha_m[beta] + 1);
		for (int kkk = 1; kkk < bicore_index_v[beta].size(); kkk++) {
			bicore_index_v[beta][kkk] = new bicore_index_block;
		}
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		for (int beta = 1; beta < g.right_index[v].size(); beta++) {
			int alpha = g.right_index[v][beta];
			bicore_index_v[beta][alpha]->nodeset.push_back(v);
		}
	}
	for (int beta = 1; beta < bicore_index_v.size(); beta++) {
		bicore_index_block* pre = NULL;
		for (int alpha = bicore_index_v[beta].size() - 1; alpha > 0; alpha--) {
			bicore_index_v[beta][alpha]->next = pre;
			if (bicore_index_v[beta][alpha]->nodeset.size() > 0) pre = bicore_index_v[beta][alpha];
		}
	}
}

void build_skyline_index(BiGraph&g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, vector<skyline_block*>& skyline_index_u, vector<skyline_block*>& skyline_index_v){
    skyline_index_u.resize(g.num_v1 + 1);
    skyline_index_v.resize(g.num_v2 + 1);
    for (vid_t w = 0; w < g.num_v1 + 1; w++){
        skyline_index_u[w] = new skyline_block;
    }
    for (vid_t w = 0; w < g.num_v2 + 1; w++){
        skyline_index_v[w] = new skyline_block;
    }
    // build right index
    for (int alpha = 1; alpha < bicore_index_u.size(); alpha++){
        for (int beta = bicore_index_u[alpha].size() - 1; beta > 0; beta--) {
            if (bicore_index_u[alpha][beta] != NULL) {
                for (auto i = bicore_index_u[alpha][beta]->nodeset.begin(); i != bicore_index_u[alpha][beta]->nodeset.end(); i++) {
                    if (skyline_index_u[*i]->vector_set.size() != 0){
                        if (skyline_index_u[*i]->vector_set[skyline_index_u[*i]->vector_set.size() - 2] >= alpha && skyline_index_u[*i]->vector_set[skyline_index_u[*i]->vector_set.size() - 1] >= beta) continue;
                        if (skyline_index_u[*i]->vector_set[skyline_index_u[*i]->vector_set.size() - 2] <= alpha && skyline_index_u[*i]->vector_set[skyline_index_u[*i]->vector_set.size() - 1] <= beta){
                            skyline_index_u[*i]->vector_set.pop_back();
                            skyline_index_u[*i]->vector_set.pop_back();
                        }
                    }
                    skyline_index_u[*i]->vector_set.push_back(alpha);
                    skyline_index_u[*i]->vector_set.push_back(beta);
                }
            }
        }
    }
    // build left index
	for (int beta = 1; beta < bicore_index_v.size(); beta++){
        for (int alpha = bicore_index_v[beta].size() - 1; alpha > 0; alpha--) {
            if (bicore_index_v[beta][alpha] != NULL) {
                for (auto i = bicore_index_v[beta][alpha]->nodeset.begin(); i != bicore_index_v[beta][alpha]->nodeset.end(); i++) {
                    if (skyline_index_v[*i]->vector_set.size() != 0){
                        if (skyline_index_v[*i]->vector_set[0] >= alpha && skyline_index_v[*i]->vector_set[1] >= beta) continue;
                        if (skyline_index_v[*i]->vector_set[0] <= alpha && skyline_index_v[*i]->vector_set[1] <= beta){
                            skyline_index_v[*i]->vector_set.erase(skyline_index_v[*i]->vector_set.begin());
                            skyline_index_v[*i]->vector_set.erase(skyline_index_v[*i]->vector_set.begin());
                        }
                    }
                    skyline_index_v[*i]->vector_set.insert(skyline_index_v[*i]->vector_set.begin(), beta);
                    skyline_index_v[*i]->vector_set.insert(skyline_index_v[*i]->vector_set.begin(), alpha);
                }
            }
        }
    }
}

// arrange the node based on the skyline point
void build_bicore_index_space_saver(BiGraph&g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v) {
	bicore_index_u.clear(); bicore_index_u.resize(g.v1_max_degree + 1); bicore_index_v.clear(); bicore_index_v.resize(g.v2_max_degree + 1);
	// build left
	vector<int> beta_m; beta_m.resize(g.v1_max_degree + 1);
	for (vid_t u = 0; u < g.num_v1; u++) {
		for (int alpha = 1; alpha < g.left_index[u].size(); alpha++) {
			int beta = g.left_index[u][alpha];
			if (beta_m[alpha] < beta) beta_m[alpha] = beta;
		}
	}
	for (int alpha = 1; alpha < beta_m.size(); alpha++) {
		bicore_index_u[alpha].resize(beta_m[alpha] + 1);
		for (int kkk = 1; kkk < bicore_index_u[alpha].size(); kkk++) {
			bicore_index_u[alpha][kkk] = new bicore_index_block;
		}
	}
	for (vid_t u = 0; u < g.num_v1; u++) {
		int alpha = g.left_index[u].size() - 1;
		if (alpha <= 0) continue;
		int pre = g.left_index[u][alpha];
		bicore_index_u[alpha][pre]->nodeset.push_back(u);
		alpha--;
		while (alpha > 0) {
			int beta = g.left_index[u][alpha];
			if (beta > pre) {
				pre = beta;
				bicore_index_u[alpha][pre]->nodeset.push_back(u);
			}
			alpha--;
		}
	}
	for (int alpha = 1; alpha < bicore_index_u.size(); alpha++) {
		bicore_index_block* pre = NULL;
		for (int beta = bicore_index_u[alpha].size() - 1; beta > 0; beta--) {
			bicore_index_u[alpha][beta]->next = pre;
			if (bicore_index_u[alpha][beta]->nodeset.size() > 0) pre = bicore_index_u[alpha][beta];
		}
	}
	//build right
	vector<int> alpha_m; alpha_m.resize(g.v2_max_degree + 1);
	for (vid_t v = 0; v < g.num_v2; v++) {
		for (int beta = 1; beta < g.right_index[v].size(); beta++) {
			int alpha = g.right_index[v][beta];
			if (alpha_m[beta] < alpha) alpha_m[beta] = alpha;
		}
	}
	for (int beta = 1; beta < alpha_m.size(); beta++) {
		bicore_index_v[beta].resize(alpha_m[beta] + 1);
		for (int kkk = 1; kkk < bicore_index_v[beta].size(); kkk++) {
			bicore_index_v[beta][kkk] = new bicore_index_block;
		}
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		int beta = g.right_index[v].size() - 1;
		if (beta <= 0) continue;
		int pre = g.right_index[v][beta];
		bicore_index_v[beta][pre]->nodeset.push_back(v);
		beta--;
		while (beta > 0) {
			int alpha = g.right_index[v][beta];
			if (alpha > pre) {
				pre = alpha;
				bicore_index_v[beta][pre]->nodeset.push_back(v);
			}
			beta--;
		}
	}
	for (int beta = 1; beta < bicore_index_v.size(); beta++) {
		bicore_index_block* pre = NULL;
		for (int alpha = bicore_index_v[beta].size() - 1; alpha > 0; alpha--) {
			bicore_index_v[beta][alpha]->next = pre;
			if (bicore_index_v[beta][alpha]->nodeset.size() > 0) pre = bicore_index_v[beta][alpha];
		}
	}
}


bool isNodesetEqual(unordered_set<vid_t>& set1, unordered_set<vid_t>& set2) {
	if (set1.size() != set2.size()) {
		return false;
	}
	for (auto it = set1.begin(); it != set1.end(); it++) {
		auto got = set2.find(*it);
		if (got == set2.end()) {
			return false;
		}
	}
	return true;
}

void dynamic_construct_dual_pointer(vector<vector<bicore_index_block_dual_pointer*>>& bicore_index) {
	vector<pair<int, unordered_set<vid_t>>> alpha_list; alpha_list.resize(bicore_index.size());
	for (int i = alpha_list.size() - 1; i > 0; i--) {
		alpha_list[i].first = bicore_index[i].size() - 1;
		alpha_list[i].second = (bicore_index[i][alpha_list[i].first])->nodeset;
		if (i == alpha_list.size() - 1) continue;
		if (alpha_list[i + 1].first == alpha_list[i].first) {
			alpha_list[i].second.insert(alpha_list[i + 1].second.begin(), alpha_list[i + 1].second.end());
			bicore_index[i][alpha_list[i].first]->vertical_pointer = bicore_index[i + 1][alpha_list[i + 1].first];
		}
	}
	int alpha_list_size = alpha_list.size();
	while (true) {
		//////////////////
		if (alpha_list[1].first == 25) {
			int opt = 0;
		}
		//////////////////
		// stop condition
		if (alpha_list[1].first == 1) {
			break;
		}
		// set sentinel
		alpha_list[0].first = alpha_list[1].first;
		for (int alpha = 1; alpha < alpha_list_size; alpha++) {
			// beta == 1
			if (alpha_list[alpha].first == 1) {
				// nothing to do in this loop
				// delete useless node set to save space and break
				if (alpha_list[alpha - 1].first == 1) {
					alpha_list[alpha].second.clear();
					break;
				}
			}
			// trigger arrived
			else if (alpha_list[alpha].first == alpha_list[alpha - 1].first) {
				// boundary condition: expand nodeset
				if (alpha + 1 == alpha_list_size || alpha_list[alpha + 1].first < alpha_list[alpha].first - 1) {
					bicore_index_block_dual_pointer* block = bicore_index[alpha][alpha_list[alpha].first - 1];
					alpha_list[alpha].second.insert(block->nodeset.begin(), block->nodeset.end());
					alpha_list[alpha].first--;
				}
				// cross point
				else if (alpha_list[alpha + 1].first == alpha_list[alpha].first - 1) {
					// set vertical_pointer
					if (!isNodesetEqual(alpha_list[alpha + 1].second, alpha_list[alpha].second)) {
						bicore_index[alpha][alpha_list[alpha].first - 1]->vertical_pointer = bicore_index[alpha + 1][alpha_list[alpha].first - 1];
					}
					bicore_index_block_dual_pointer* block = bicore_index[alpha][alpha_list[alpha].first - 1];
					alpha_list[alpha].second.insert(block->nodeset.begin(), block->nodeset.end());
					block = bicore_index[alpha + 1][alpha_list[alpha].first - 1];
					alpha_list[alpha].second.insert(block->nodeset.begin(), block->nodeset.end());
					alpha_list[alpha].first--;
				}
			}
			// hold on waiting for trigger
			else {
				// do nothing
			}
		}
	}
}

void build_bicore_index_space_saver_dual_pointer(BiGraph&g, vector<vector<bicore_index_block_dual_pointer*>>& bicore_index_u, vector<vector<bicore_index_block_dual_pointer*>>& bicore_index_v) {
	bicore_index_u.clear(); bicore_index_u.resize(g.v1_max_degree + 1); bicore_index_v.clear(); bicore_index_v.resize(g.v2_max_degree + 1);
	// build left
	vector<int> beta_m; beta_m.resize(g.v1_max_degree + 1);
	for (vid_t u = 0; u < g.num_v1; u++) {
		for (int alpha = 1; alpha < g.left_index[u].size(); alpha++) {
			int beta = g.left_index[u][alpha];
			if (beta_m[alpha] < beta) beta_m[alpha] = beta;
		}
	}
	for (int alpha = 1; alpha < beta_m.size(); alpha++) {
		bicore_index_u[alpha].resize(beta_m[alpha] + 1);
		for (int kkk = 1; kkk < bicore_index_u[alpha].size(); kkk++) {
			bicore_index_u[alpha][kkk] = new bicore_index_block_dual_pointer;
		}
	}
	for (vid_t u = 0; u < g.num_v1; u++) {
		int alpha = g.left_index[u].size() - 1;
		if (alpha <= 0) continue;
		int pre = g.left_index[u][alpha];
		bicore_index_u[alpha][pre]->nodeset.insert(u);
		alpha--;
		while (alpha > 0) {
			int beta = g.left_index[u][alpha];
			if (beta > pre) {
				pre = beta;
				bicore_index_u[alpha][pre]->nodeset.insert(u);
			}
			alpha--;
		}
	}
	for (int alpha = 1; alpha < bicore_index_u.size(); alpha++) {
		bicore_index_block_dual_pointer* pre = NULL;
		for (int beta = bicore_index_u[alpha].size() - 1; beta > 0; beta--) {
			bicore_index_u[alpha][beta]->horizontal_pointer = pre;
			pre = bicore_index_u[alpha][beta];
		}
	}
	//build right
	vector<int> alpha_m; alpha_m.resize(g.v2_max_degree + 1);
	for (vid_t v = 0; v < g.num_v2; v++) {
		for (int beta = 1; beta < g.right_index[v].size(); beta++) {
			int alpha = g.right_index[v][beta];
			if (alpha_m[beta] < alpha) alpha_m[beta] = alpha;
		}
	}
	for (int beta = 1; beta < alpha_m.size(); beta++) {
		bicore_index_v[beta].resize(alpha_m[beta] + 1);
		for (int kkk = 1; kkk < bicore_index_v[beta].size(); kkk++) {
			bicore_index_v[beta][kkk] = new bicore_index_block_dual_pointer;
		}
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		int beta = g.right_index[v].size() - 1;
		if (beta <= 0) continue;
		int pre = g.right_index[v][beta];
		bicore_index_v[beta][pre]->nodeset.insert(v);
		beta--;
		while (beta > 0) {
			int alpha = g.right_index[v][beta];
			if (alpha > pre) {
				pre = alpha;
				bicore_index_v[beta][pre]->nodeset.insert(v);
			}
			beta--;
		}
	}
	for (int beta = 1; beta < bicore_index_v.size(); beta++) {
		bicore_index_block_dual_pointer* pre = NULL;
		for (int alpha = bicore_index_v[beta].size() - 1; alpha > 0; alpha--) {
			bicore_index_v[beta][alpha]->horizontal_pointer = pre;
			pre = bicore_index_v[beta][alpha];
		}
	}
	// dynamic construct
	dynamic_construct_dual_pointer(bicore_index_u);
	dynamic_construct_dual_pointer(bicore_index_v);
}

void retrieve_via_bicore_index_space_saver_dual_pointer(BiGraph& g, vector<vector<bicore_index_block_dual_pointer*>>& bicore_index_u, vector<vector<bicore_index_block_dual_pointer*>>& bicore_index_v,
	vector<bool>& left_node, vector<bool>& right_node, int alpha, int beta) {
	left_node.clear(); right_node.clear();
	left_node.resize(g.num_v1); right_node.resize(g.num_v2);
	fill_n(left_node.begin(), left_node.size(), false);
	fill_n(right_node.begin(), right_node.size(), false);
	if (bicore_index_u.size() <= alpha) return;
	if (bicore_index_u[alpha].size() <= beta) return;
	vector<bicore_index_block_dual_pointer*> queue;
	bicore_index_block_dual_pointer* block = bicore_index_u[alpha][beta];
	while (block != NULL) {
		if (block->vertical_pointer != NULL) {
			queue.push_back(block->vertical_pointer);
		}
		for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
			left_node[*i] = true;
		}
		block = block->horizontal_pointer;
	}
	while (!queue.empty()) {
		block = queue.back();
		queue.pop_back();
		while (block != NULL) {
			for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
				left_node[*i] = true;
			}
			block = block->vertical_pointer;
		}
	}
	block = bicore_index_v[beta][alpha];
	while (block != NULL) {
		if (block->vertical_pointer != NULL) {
			queue.push_back(block->vertical_pointer);
		}
		for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
			right_node[*i] = true;
		}
		block = block->horizontal_pointer;
	}
	while (!queue.empty()) {
		block = queue.back();
		queue.pop_back();
		while (block != NULL) {
			for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
				right_node[*i] = true;
			}
			block = block->vertical_pointer;
		}
	}
}

void retrieve_via_bicore_index_space_saver(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v,
	vector<bool>& left_node, vector<bool>& right_node, int alpha, int beta) {
	left_node.clear(); right_node.clear();
	left_node.resize(g.num_v1); right_node.resize(g.num_v2);
	fill_n(left_node.begin(), left_node.size(), false);
	fill_n(right_node.begin(), right_node.size(), false);
	if (bicore_index_u.size() <= alpha) return;
	if (bicore_index_u[alpha].size() <= beta) return;
	int alpha_ = alpha;
	int beta_ = beta;
	while (bicore_index_u.size() > alpha_) {
		if (bicore_index_u[alpha_].size() <= beta_) {
			break;
		}
		bicore_index_block* block = bicore_index_u[alpha_][beta_];
		while (block != NULL) {
			for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
				left_node[*i] = true;
			}
			block = block->next;
		}
		alpha_++;
	}
	alpha_ = alpha;
	beta_ = beta;
	while (bicore_index_v.size() > beta_) {
		if (bicore_index_v[beta_].size() <= alpha_) {
			break;
		}
		bicore_index_block* block = bicore_index_v[beta_][alpha_];
		while (block != NULL) {
			for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
				right_node[*i] = true;
			}
			block = block->next;
		}
		beta_++;
	}
}

void retrieve_via_bicore_index(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v,
	vector<bool>& left_node, vector<bool>& right_node, int alpha, int beta) {
	//left_node.clear(); right_node.clear();
	//left_node.resize(g.num_v1); right_node.resize(g.num_v2);
	fill_n(left_node.begin(), left_node.size(), false);
	fill_n(right_node.begin(), right_node.size(), false);
	if (bicore_index_u.size() <= alpha) return;
	if (bicore_index_u[alpha].size() <= beta) return;
	bicore_index_block* block = bicore_index_u[alpha][beta];
	while (block != NULL) {
		for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
			left_node[*i] = true;
		}
		block = block->next;
	}
	block = bicore_index_v[beta][alpha];
	while (block != NULL) {
		for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
			right_node[*i] = true;
		}
		block = block->next;
	}
}

void retrieve_node_via_bicore_index (BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, batch_set& pair, vid_t node, bool is_left){
	if (is_left){
		for(int alpha = 1; alpha < bicore_index_u.size(); alpha++){
			bool find = false;
			for (int beta = 1; beta < bicore_index_u[alpha].size(); beta++){
				bool find_ = false;
				bicore_index_block* block = bicore_index_u[alpha][beta];
				if (block == NULL) break;
				for (auto node_ : block->nodeset){
					if(node_ == node){
						find = true;
						find_ = true;
						break;
					}
				}
				if (find && !find_){
					//for (int j = 1; j < beta; j++){
						pair.insert(make_pair(alpha , beta - 1));
					//}
					break;
				}
				if (find && block->next == NULL){
					//for (int j = 1; j <= beta; j++){
						pair.insert(make_pair(alpha , beta));
					//}
					break;
				}
			}
			if (!find) break;
		}
	}
	else{
		for (int beta = 1; beta < bicore_index_v.size(); beta++){
			bool find = false;
			for (int alpha = 1; alpha < bicore_index_v[beta].size(); alpha++){
				bool find_ = false;
				bicore_index_block* block = bicore_index_v[beta][alpha];
				if (block == NULL) break;
				for (auto node_ : block->nodeset){
					if(node_ == node){
						//for (int j = 1; j <= beta; j++){
						//	pair.insert(make_pair(alpha , j));
						//}
						find = true;
						find_ = true;
						break;
					}
				}
				if (find && !find_){
					//for (int i = 1; i < alpha; i++){
						pair.insert(make_pair(alpha - 1 , beta));
					//}
					break;
				}
				if (find && block->next == NULL){
					//for (int i = 1; i <= alpha; i++){
						pair.insert(make_pair(alpha , beta));
					//}
					break;
				}
			}
			if (!find) break;
		}
	}
}

// !important! this function will not initialize left_node and right_node
void retrieve_via_bicore_index_inverse(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v,
	vector<bool>& left_node, vector<bool>& right_node, int alpha, int beta) {
	if (bicore_index_u.size() <= alpha) return;
	if (bicore_index_u[alpha].size() <= beta) return;
	bicore_index_block* block = bicore_index_u[alpha][beta];
	while (block != NULL) {
		for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
			left_node[*i] = false;
		}
		block = block->next;
	}
	block = bicore_index_v[beta][alpha];
	while (block != NULL) {
		for (auto i = block->nodeset.begin(); i != block->nodeset.end(); i++) {
			right_node[*i] = false;
		}
		block = block->next;
	}
}

void retrieve_node_via_index_skyline (BiGraph& g, vector<skyline_block*>& skyline_index_u, vector<skyline_block*>& skyline_index_v, batch_set& pair, vid_t node, bool is_left){
	if(is_left){
		int i_pre = 0;
		skyline_block* block = skyline_index_u[node];
		for (auto i = 0; i < block->vector_set.size(); i+=2){
			int alpha = block->vector_set[i]; int beta = block->vector_set[i + 1];
			for (int alpha_ = i_pre + 1; alpha_ <= alpha; alpha_++){
			//	for (int beta = 1; beta <= j; beta++) {
					pair.insert(make_pair(alpha_, beta));
			//	}
			}
			i_pre = alpha;
		}
	}
	else{
		int i_pre = 0;
		skyline_block* block = skyline_index_v[node];
		for (auto i = 0; i < block->vector_set.size(); i+=2){
			int alpha = block->vector_set[i]; int beta = block->vector_set[i + 1];
			for (int alpha_ = i_pre + 1; alpha_ <= alpha; alpha_++){
			//	for (int beta = 1; beta <= j; beta++) {
					pair.insert(make_pair(alpha_, beta));
			//	}
			}
			i_pre = alpha;
		}
	}
}

int binseaech(skyline_block* array, int len, int target)
{
	//for (auto i : array->vector_set) cout << i << endl;
	int left = 0, right = len/2 - 1;
	int mid;
	while (left<=right)
	{
		mid = (left + right) / 2;
		if (target == array->vector_set[2*mid])
		{
			return mid;
		}
		else if (target > array->vector_set[2*mid])
		{
			left = mid + 1;
		}
		else if (target < array->vector_set[2*mid])
		{
			right = mid - 1;
		}
	}
	return left;
}


void retrieve_bicore_via_bicore_number (BiGraph& g, vector<skyline_block*>& skyline_index_u, vector<skyline_block*>& skyline_index_v,
	vector<bool>& left_node, vector<bool>& right_node, int alpha, int beta){
	fill_n(left_node.begin(), left_node.size(), false);
	fill_n(right_node.begin(), right_node.size(), false);
	for (int u = 0; u < skyline_index_u.size(); u++){
		skyline_block* block = skyline_index_u[u];
		if (block->vector_set.size() == 0) continue;
		if (block->vector_set[block->vector_set.size() - 2] < alpha || block->vector_set[1] < beta) continue;
		int result = binseaech(block, block->vector_set.size(), alpha);
		if (block->vector_set[result*2 + 1] >= beta){
			left_node[u] = true;
		}
	}
	for (int v = 0; v < skyline_index_v.size(); v++){
		skyline_block* block = skyline_index_v[v];
		if (block->vector_set.size() == 0) continue;
		if (block->vector_set[block->vector_set.size() - 2] < alpha || block->vector_set[1] < beta) continue;
		int result = binseaech(block, block->vector_set.size(), alpha);
		if (block->vector_set[result*2 + 1] >= beta){
			right_node[v] = true;
		}
	}
}

void retrieve_hierarchy_via_index_skyline (BiGraph& g, vector<skyline_block*>& skyline_index_u, vector<skyline_block*>& skyline_index_v, 
	vector<int>& left_node, vector<int>& right_node, queue<int> node_queue, queue<bool> is_left_queue, int alpha, int beta){
	while(is_left_queue.size() > 0){
		vid_t node = node_queue.front(); node_queue.pop();
		bool is_left = is_left_queue.front(); is_left_queue.pop();
		if(is_left){
			skyline_block* block = skyline_index_u[node];
			int size = block->vector_set.size();
			if (size != 0){
				if (block->vector_set[size - 2] >= alpha && block->vector_set[1] >= beta){
					int result = binseaech(block, size, alpha);
					if (block->vector_set[result*2 + 1] >= beta){
						left_node[node] = 1;
						for (vid_t v : g.neighbor_v1[node]){
							if (right_node[v] == 0){
								right_node[v] = 2;
								node_queue.push(v);
								is_left_queue.push(false);
							}
						}
					}
				}
			}
		}
		else{
			skyline_block* block = skyline_index_v[node];
			int size = block->vector_set.size();
			if (size != 0){
				if (block->vector_set[size - 2] >= alpha && block->vector_set[1] >= beta){
					int result = binseaech(block, size, alpha);
					if (block->vector_set[result*2 + 1] >= beta){
						right_node[node] = 1;
						for (vid_t u : g.neighbor_v2[node]){
							if (left_node[u] == 0){
								left_node[u] = 2;
								node_queue.push(u);
								is_left_queue.push(true);
							}
						}
					}
				}
			}
		}
	}
}

void hierarchy_bfs(BiGraph& g, vector<bool>& left_node, vector<bool>& right_node, set<int>& left_set, set<int>& right_set, vector<int>& left_visit, vector<int>& right_visit, queue<int>& node_queue, queue<bool>& is_left_queue){
	vid_t node = node_queue.front(); bool is_left = is_left_queue.front();
	node_queue.pop(); is_left_queue.pop();
	if(is_left){
		if (left_node[node]){
			left_set.insert(node);
			for(vid_t v : g.neighbor_v1[node]){
				if(right_visit[v] == 0) {
					right_visit[v] = 1;
					node_queue.push(v); is_left_queue.push(false);
				}
			}
		}
	}
	else{
		if (right_node[node]){
			right_set.insert(node);
			for(vid_t u : g.neighbor_v2[node]){
				if (left_visit[u] == 0) {
					left_visit[u] = 1;
					node_queue.push(u); is_left_queue.push(true);
				}
			}
		}
	}
}

void retrieve_hierarchy_via_bicore_index(BiGraph& g, vector<vector<bicore_index_block*>>& bicore_index_u, vector<vector<bicore_index_block*>>& bicore_index_v, 
	vector<bool>& left_node, vector<bool>& right_node, set<int>& left_set, set<int>& right_set, vector<int>& left_visit, vector<int>& right_visit, vid_t node, bool is_left, int alpha, int beta){
		retrieve_via_bicore_index(g, bicore_index_u, bicore_index_v, left_node, right_node, alpha, beta);
		queue<int> node_queue; queue<bool> is_left_queue;
		node_queue.push(node); is_left_queue.push(is_left);
		left_visit[node] = 1;
		while(node_queue.size() > 0){
			hierarchy_bfs(g, left_node, right_node, left_set, right_set, left_visit, right_visit, node_queue, is_left_queue);
		}
	}

void inv(BiGraph& g) {
	swap(g.degree_v1, g.degree_v2);
	swap(g.left_index, g.right_index);
	swap(g.num_v1, g.num_v2);
	swap(g.neighbor_v1, g.neighbor_v2);
	swap(g.neighborHash_v1, g.neighborHash_v2);
	swap(g.v1_max_degree, g.v2_max_degree);
	swap(g.left_delete, g.right_delete);
}
