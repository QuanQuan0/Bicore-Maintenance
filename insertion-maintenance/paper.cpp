#include "paper.h"
#include <algorithm>
#include <ctime>
#include <iostream>
#include <set>
#include <cmath>
#include <list>

using namespace std;

bool DELETEC = false;
bool APC = false;
long long DELETECOUNT = 0;
long long APCOUNT = 0;

void print_index_detail(BiGraph& g) {
	for (vid_t u = 0; u < g.num_v1; u++) {
		cout << "u" << u + 1 << " : ";
		for (int i = 1; i < g.left_index[u].size(); i++) {
			cout << g.left_index[u][i] << " ";
		}
		cout << endl;
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		cout << "v" << v + 1 << " : ";
		for (int i = 1; i < g.right_index[v].size(); i++) {
			cout << g.right_index[v][i] << " ";
		}
		cout << endl;
	}
}

void binsort_paper(vector<list<vid_t>>& left_bin, vector<list<vid_t>>& right_bin,
	vector<list<vid_t>::iterator>& left_list_pointer, vector<list<vid_t>::iterator>& right_list_pointer, BiGraph& tmp) {
	left_bin.clear(); right_bin.clear(); left_list_pointer.clear(); right_list_pointer.clear();
	left_bin.resize(tmp.v1_max_degree + 1);
	left_list_pointer.resize(tmp.getV1Num());
	for (int i = 0; i < tmp.getV1Num(); i++) {
		left_bin[tmp.getV1Degree(i)].push_back(i);
		left_list_pointer[i] = left_bin[tmp.getV1Degree(i)].end();
		left_list_pointer[i]--;
	}
	right_list_pointer.resize(tmp.getV2Num());
	right_bin.resize(tmp.v2_max_degree + 1);
	for (int i = 0; i < tmp.getV2Num(); i++) {
		right_bin[tmp.getV2Degree(i)].push_back(i);
		right_list_pointer[i] = right_bin[tmp.getV2Degree(i)].end();
		right_list_pointer[i]--;
	}
}

void reorder_left_paper(vector<list<vid_t>>& bin, vector<list<vid_t>::iterator>& list_pointer, int new_degree, int old_degree, vid_t vertex) {
	bin[old_degree].erase(list_pointer[vertex]);
	bin[new_degree].push_back(vertex);
	list_pointer[vertex] = bin[new_degree].end();
	list_pointer[vertex]--;
}

void reorder_right_paper(vector<list<vid_t>>& bin, vector<list<vid_t>::iterator>& list_pointer, int new_degree, int old_degree, vid_t vertex) {
	bin[old_degree].erase(list_pointer[vertex]);
	bin[new_degree].push_back(vertex);
	list_pointer[vertex] = bin[new_degree].end();
	list_pointer[vertex]--;
}

bool compare_deg_array(const pair<int, int>& p1, const pair<int, int>& p2) {
	return p1.second < p2.second;
}

int maxkcorenumber(BiGraph& g) {
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;

	// definition of left dyn array
	vector<int> left_k_x_index; left_k_x_index.resize(g.v1_max_degree + 2);
	fill_n(left_k_x_index.begin(), left_k_x_index.size(), -1); left_k_x_index[0] = 0; left_k_x_index[g.v1_max_degree + 1] = g.num_v1;
	vector<int> left_vertices_index; left_vertices_index.resize(g.num_v1);
	vector<pair<int, int>> left_degree_array; left_degree_array.resize(g.num_v1);
	// definition of right dyn array
	vector<int> right_k_x_index; right_k_x_index.resize(g.v2_max_degree + 2);
	fill_n(right_k_x_index.begin(), right_k_x_index.size(), -1); right_k_x_index[0] = 0; right_k_x_index[g.v2_max_degree + 1] = g.num_v2;
	vector<int> right_vertices_index; right_vertices_index.resize(g.num_v2);
	vector<pair<int, int>> degree_array; degree_array.resize(g.num_v2);
	// end of definition

	// initialize left dyn array
	for (vid_t u = 0; u < g.num_v1; u++) {
		left_degree_array[u] = make_pair(u, g.degree_v1[u]);
	}
	sort(left_degree_array.begin(), left_degree_array.end(), compare_deg_array);
	for (int i = 0; i < left_degree_array.size(); i++) {
		if (left_k_x_index[left_degree_array[i].second] == -1) {
			left_k_x_index[left_degree_array[i].second] = i;
		}
		left_vertices_index[left_degree_array[i].first] = i;
	}
	for (int i = left_k_x_index.size() - 1; i >= 0; i--) {
		if (left_k_x_index[i] == -1) {
			left_k_x_index[i] = left_k_x_index[i + 1];
		}
	}
	// initialize right dyn array
	for (vid_t v = 0; v < g.num_v2; v++) {
		degree_array[v] = make_pair(v, g.degree_v2[v]);
	}
	sort(degree_array.begin(), degree_array.end(), compare_deg_array);
	for (int i = 0; i < degree_array.size(); i++) {
		if (right_k_x_index[degree_array[i].second] == -1) {
			right_k_x_index[degree_array[i].second] = i;
		}
		right_vertices_index[degree_array[i].first] = i;
	}
	for (int i = right_k_x_index.size() - 1; i >= 0; i--) {
		if (right_k_x_index[i] == -1) {
			right_k_x_index[i] = right_k_x_index[i + 1];
		}
	}
	// end of initialization
	int kc = 1;
	for (kc = 1; kc <= g.v2_max_degree + 1; kc++) {
		if (left_k_x_index[kc - 1] == g.num_v1) break;
		if (right_k_x_index[kc - 1] == g.num_v2) break;
		for (int i = left_k_x_index[kc - 1]; i < left_k_x_index[kc]; i++) {
			if (!g.left_delete[left_degree_array[i].first]) left_vertices_to_be_peeled.push_back(left_degree_array[i].first);
		}
		for (int i = right_k_x_index[kc - 1]; i < right_k_x_index[kc]; i++) {
			if (!g.right_delete[degree_array[i].first]) right_vertices_to_be_peeled.push_back(degree_array[i].first);
		}
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			for (auto j = left_vertices_to_be_peeled.begin(); j != left_vertices_to_be_peeled.end(); j++) {
				vid_t u = *j;
				if (g.left_delete[u]) continue;
				for (int k = 0; k < g.neighbor_v1[u].size(); k++) {
					vid_t v = g.neighbor_v1[u][k];
					if (g.right_delete[v]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.degree_v2[v] == 0) {
						g.right_delete[v] = true;
					}
					if (!g.right_delete[v] && g.degree_v2[v] < kc) {
						right_vertices_to_be_peeled.push_back(v);
					}
					if (g.degree_v2[v] >= kc - 1) {
						int olddegree = g.degree_v2[v] + 1;
						vid_t stack = degree_array[right_k_x_index[olddegree]].first;
						swap(degree_array[right_k_x_index[olddegree]], degree_array[right_vertices_index[v]]);
						swap(right_vertices_index[stack], right_vertices_index[v]);
						right_k_x_index[olddegree]++;
					}
				}
				g.left_delete[u] = true;
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			for (auto j = right_vertices_to_be_peeled.begin(); j != right_vertices_to_be_peeled.end(); j++) {
				vid_t v = *j;
				if (g.right_delete[v]) continue;
				for (int k = 0; k < g.neighbor_v2[v].size(); k++) {
					vid_t u = g.neighbor_v2[v][k];
					if (g.left_delete[u]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.degree_v1[u] == 0) {
						g.left_delete[u] = true;
					}
					if (!g.left_delete[u] && g.degree_v1[u] < kc) {
						left_vertices_to_be_peeled.push_back(u);
					}
					if (g.degree_v1[u] >= kc - 1) {
						int olddegree = g.degree_v1[u] + 1;
						vid_t stack = left_degree_array[left_k_x_index[olddegree]].first;
						swap(left_degree_array[left_k_x_index[olddegree]], left_degree_array[left_vertices_index[u]]);
						swap(left_vertices_index[stack], left_vertices_index[u]);
						left_k_x_index[olddegree]++;
					}
				}
				g.right_delete[v] = true;
			}
			right_vertices_to_be_peeled.clear();
		}
	}
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}
	return kc - 2;
}

void crossUpdate(BiGraph& g, int left_k, int k_x, vid_t v) {
	for (int right_k = 1; right_k < g.right_index[v].size(); right_k++) {
		if (right_k <= k_x) {
			if (g.right_index[v][right_k] < left_k) {
				g.right_index[v][right_k] = left_k;
			}
		}
		else {
			break;
		}
	}
}

void __alphaPeel(int left_k, BiGraph& g) {
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	// definition of dyn array
	vector<int> right_k_x_index; right_k_x_index.resize(g.v2_max_degree + 2);
	fill_n(right_k_x_index.begin(), right_k_x_index.size(), -1); right_k_x_index[0] = 0; right_k_x_index[g.v2_max_degree + 1] = g.num_v2;
	vector<int> right_vertices_index; right_vertices_index.resize(g.num_v2);
	vector<pair<int, int>> degree_array; degree_array.resize(g.num_v2);
	// end of definition
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
			left_vertices_to_be_peeled.push_back(u);
		}
	}
	// initialize dyn array
	for (vid_t v = 0; v < g.num_v2; v++) {
		degree_array[v] = make_pair(v, g.degree_v2[v]);
	}
	sort(degree_array.begin(), degree_array.end(), compare_deg_array);
	for (int i = 0; i < degree_array.size(); i++) {
		if (right_k_x_index[degree_array[i].second] == -1) {
			right_k_x_index[degree_array[i].second] = i;
		}
		right_vertices_index[degree_array[i].first] = i;
	}
	for (int i = right_k_x_index.size() - 1; i >= 0; i--) {
		if (right_k_x_index[i] == -1) {
			right_k_x_index[i] = right_k_x_index[i + 1];
		}
	}
	// end of initialization
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		if (right_k_x_index[right_k - 1] == g.num_v2) break;
		for (int i = right_k_x_index[right_k - 1]; i < right_k_x_index[right_k]; i++) {
			right_vertices_to_be_peeled.push_back(degree_array[i].first);
		}
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			for (auto j = left_vertices_to_be_peeled.begin(); j != left_vertices_to_be_peeled.end(); j++) {
				vid_t u = *j;
				if (g.left_delete[u]) continue;
				for (int k = 0; k < g.neighbor_v1[u].size(); k++) {
					vid_t v = g.neighbor_v1[u][k];
					if (g.right_delete[v]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV2Degree(v) == 0 && right_k - 1 > 0) {
						g.right_delete[v] = true;
					}
					if (g.degree_v2[v] < right_k && !g.right_delete[v]) {
						right_vertices_to_be_peeled.push_back(v);
					}
					if (g.degree_v2[v] >= right_k - 1) {
						int olddegree = g.degree_v2[v] + 1;
						vid_t stack = degree_array[right_k_x_index[olddegree]].first;
						swap(degree_array[right_k_x_index[olddegree]], degree_array[right_vertices_index[v]]);
						swap(right_vertices_index[stack], right_vertices_index[v]);
						right_k_x_index[olddegree]++;
					}
				}
				g.left_delete[u] = true;
				if (right_k - 1 > 0) {
					g.left_index[u][left_k] = right_k - 1;
				}
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			for (auto j = right_vertices_to_be_peeled.begin(); j != right_vertices_to_be_peeled.end(); j++) {
				vid_t v = *j;
				if (g.right_delete[v]) continue;
				for (int k = 0; k < g.neighbor_v2[v].size(); k++) {
					vid_t u = g.neighbor_v2[v][k];
					if (g.left_delete[u]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV1Degree(u) == 0 && right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
						g.left_delete[u] = true;
					}
					if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.right_delete[v] = true;
			}
			right_vertices_to_be_peeled.clear();
		}
	}
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}
}

void alphaPeel(int left_k, BiGraph& tmp, BiGraph& g) {
	BiGraph newtmp;
	vector<list<vid_t>> left_bin;
	vector<list<vid_t>::iterator> left_list_pointer;
	vector<list<vid_t>> right_bin;
	vector<list<vid_t>::iterator> right_list_pointer;
	binsort_paper(left_bin, right_bin, left_list_pointer, right_list_pointer, tmp);
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		bool left_unchange = false;
		bool right_unchange = false;
		while (!left_unchange || !right_unchange) {
			left_unchange = true;
			right_unchange = true;
			// peel left
			for (int i = 1; i < left_k; i++) {
				for (auto j = left_bin[i].begin(); j != left_bin[i].end();) {
					vid_t u = *j;
					vector<vid_t> neigh = tmp.getV1Neighbors(u);
					int old_degree = tmp.getV1Degree(u); // should equal to i
					for (int k = 0; k < old_degree; k++) {
						vid_t v = neigh[k];
						int old_deg = tmp.getV2Degree(v); // should equal to tmp.getV2Degree(v) + 1
						tmp.deleteEdge(u, v);
						reorder_right_paper(right_bin, right_list_pointer, tmp.getV2Degree(v), old_deg, v);
					}
					++j; // j becomes unavailable after reordering
					reorder_left_paper(left_bin, left_list_pointer, tmp.getV1Degree(u), old_degree, u);
					if (right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
					}
					left_unchange = false;
				}
			}
			// peel right
			for (int i = 1; i < right_k; i++) {
				for (auto j = right_bin[i].begin(); j != right_bin[i].end(); ) {
					vid_t v = *j;
					vector<vid_t> neigh = tmp.getV2Neighbors(v);
					int old_degree = tmp.getV2Degree(v);
					for (int k = 0; k < old_degree; k++) {
						vid_t u = neigh[k];
						int old_deg = tmp.getV1Degree(u);
						tmp.deleteEdge(u, v);
						reorder_left_paper(left_bin, left_list_pointer, tmp.getV1Degree(u), old_deg, u);
						// degree = 0 will not be considered in next round
						if (tmp.getV1Degree(u) == 0 && right_k - 1 > 0) {
							g.left_index[u][left_k] = right_k - 1;
						}
					}
					++j;
					reorder_right_paper(right_bin, right_list_pointer, tmp.getV2Degree(v), old_degree, v);
					right_unchange = false;
				}
			}
		}
		if (right_k == 1) {
			newtmp = tmp;
		}
	}
	tmp = newtmp;
}

void __alphaCopyPeel(int left_k, BiGraph& g) {
	vector<bool> left_deletion_next_round;
	vector<bool> right_deletion_next_round;
	vector<int> left_degree_next_round;
	vector<int> right_degree_next_round;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	// definition of dyn array
	vector<int> right_k_x_index; right_k_x_index.resize(g.v2_max_degree + 2);
	fill_n(right_k_x_index.begin(), right_k_x_index.size(), -1); right_k_x_index[0] = 0; right_k_x_index[g.v2_max_degree + 1] = g.num_v2;
	vector<int> right_vertices_index; right_vertices_index.resize(g.num_v2);
	vector<pair<int, int>> degree_array; degree_array.resize(g.num_v2);
	// end of definition
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
			left_vertices_to_be_peeled.push_back(u);
		}
	}
	// initialize dyn array
	for (vid_t v = 0; v < g.num_v2; v++) {
		degree_array[v] = make_pair(v, g.degree_v2[v]);
	}
	sort(degree_array.begin(), degree_array.end(), compare_deg_array);
	for (int i = 0; i < degree_array.size(); i++) {
		if (right_k_x_index[degree_array[i].second] == -1) {
			right_k_x_index[degree_array[i].second] = i;
		}
		right_vertices_index[degree_array[i].first] = i;
	}
	for (int i = right_k_x_index.size() - 1; i >= 0; i--) {
		if (right_k_x_index[i] == -1) {
			right_k_x_index[i] = right_k_x_index[i + 1];
		}
	}
	// end of initialization
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		if (right_k_x_index[right_k - 1] == g.num_v2) break;
		for (int i = right_k_x_index[right_k - 1]; i < right_k_x_index[right_k]; i++) {
			right_vertices_to_be_peeled.push_back(degree_array[i].first);
		}
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			for (auto j = left_vertices_to_be_peeled.begin(); j != left_vertices_to_be_peeled.end(); j++) {
				vid_t u = *j;
				if (g.left_delete[u]) continue;
				for (int k = 0; k < g.neighbor_v1[u].size(); k++) {
					vid_t v = g.neighbor_v1[u][k];
					if (g.right_delete[v]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV2Degree(v) == 0 && right_k - 1 > 0) {
						crossUpdate(g, left_k, right_k - 1, v);
						g.right_delete[v] = true;
					}
					if (g.degree_v2[v] < right_k && !g.right_delete[v]) {
						right_vertices_to_be_peeled.push_back(v);
					}
					if (g.degree_v2[v] >= right_k - 1) {
						int olddegree = g.degree_v2[v] + 1;
						vid_t stack = degree_array[right_k_x_index[olddegree]].first;
						swap(degree_array[right_k_x_index[olddegree]], degree_array[right_vertices_index[v]]);
						swap(right_vertices_index[stack], right_vertices_index[v]);
						right_k_x_index[olddegree]++;
					}
				}
				g.left_delete[u] = true;
				if (right_k - 1 > 0) {
					g.left_index[u][left_k] = right_k - 1;
				}
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			for (auto j = right_vertices_to_be_peeled.begin(); j != right_vertices_to_be_peeled.end(); j++) {
				vid_t v = *j;
				if (g.right_delete[v]) continue;
				for (int k = 0; k < g.neighbor_v2[v].size(); k++) {
					vid_t u = g.neighbor_v2[v][k];
					if (g.left_delete[u]) continue;
					g.degree_v2[v]--;
					g.degree_v1[u]--;
					if (g.getV1Degree(u) == 0 && right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
						g.left_delete[u] = true;
					}
					if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.right_delete[v] = true;
				if (right_k - 1 > 0) {
					crossUpdate(g, left_k, right_k - 1, v);
				}
			}
			right_vertices_to_be_peeled.clear();
		}
		if (right_k == 1) {
			left_degree_next_round = g.degree_v1;
			right_degree_next_round = g.degree_v2;
			left_deletion_next_round = g.left_delete;
			right_deletion_next_round = g.right_delete;
		}
	}
	g.degree_v1 = left_degree_next_round;
	g.degree_v2 = right_degree_next_round;
	g.left_delete = left_deletion_next_round;
	g.right_delete = right_deletion_next_round;
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}

void alphaCopyPeel(int left_k, BiGraph& tmp, BiGraph& g) {
	BiGraph newtmp;
	vector<list<vid_t>> left_bin;
	vector<list<vid_t>::iterator> left_list_pointer;
	vector<list<vid_t>> right_bin;
	vector<list<vid_t>::iterator> right_list_pointer;
	binsort_paper(left_bin, right_bin, left_list_pointer, right_list_pointer, tmp);
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		bool left_unchange = false;
		bool right_unchange = false;
		while (!left_unchange || !right_unchange) {
			left_unchange = true;
			right_unchange = true;
			// peel left
			for (int i = 1; i < left_k; i++) {
				for (auto j = left_bin[i].begin(); j != left_bin[i].end();) {
					vid_t u = *j;
					vector<vid_t> neigh = tmp.getV1Neighbors(u);
					int old_degree = tmp.getV1Degree(u); // should equal to i
					for (int k = 0; k < old_degree; k++) {
						vid_t v = neigh[k];
						int old_deg = tmp.getV2Degree(v); // should equal to tmp.getV2Degree(v) + 1
						tmp.deleteEdge(u, v);
						reorder_right_paper(right_bin, right_list_pointer, tmp.getV2Degree(v), old_deg, v);
						// degree = 0 will not be considered in next round
						if (tmp.getV2Degree(v) == 0 && right_k - 1 > 0) {
							crossUpdate(g, left_k, right_k - 1, v);
						}
					}
					++j; // j becomes unavailable after reordering
					reorder_left_paper(left_bin, left_list_pointer, tmp.getV1Degree(u), old_degree, u);
					if (right_k - 1 > 0) {
						g.left_index[u][left_k] = right_k - 1;
					}
					left_unchange = false;
				}
			}
			// peel right
			for (int i = 1; i < right_k; i++) {
				for (auto j = right_bin[i].begin(); j != right_bin[i].end(); ) {
					vid_t v = *j;
					vector<vid_t> neigh = tmp.getV2Neighbors(v);
					int old_degree = tmp.getV2Degree(v);
					for (int k = 0; k < old_degree; k++) {
						vid_t u = neigh[k];
						int old_deg = tmp.getV1Degree(u);
						tmp.deleteEdge(u, v);
						reorder_left_paper(left_bin, left_list_pointer, tmp.getV1Degree(u), old_deg, u);
						// degree = 0 will not be considered in next round
						if (tmp.getV1Degree(u) == 0 && right_k - 1 > 0) {
							g.left_index[u][left_k] = right_k - 1;
						}
					}
					++j;
					reorder_right_paper(right_bin, right_list_pointer, tmp.getV2Degree(v), old_degree, v);
					if (right_k - 1 > 0) {
						crossUpdate(g, left_k, right_k - 1, v);
					}
					right_unchange = false;
				}
			}
		}
		if (right_k == 1) {
			newtmp = tmp;
		}
	}
	tmp = newtmp;
}

void crossUpdate_for_kcore(BiGraph& g, int alpha, int k_x, vid_t v) {
	for (int beta = k_x; beta > 0; beta--) {
		if (g.right_index[v][beta] < alpha) {
			g.right_index[v][beta] = alpha;
		}
		else {
			break;
		}
	}
}

void alphaCopyPeel_for_kcore(int left_k, BiGraph& g) {
	int dd_;
	int pre_left_k_ = left_k - 1;
	vector<bool> left_deletion_next_round;
	vector<bool> right_deletion_next_round;
	vector<int> left_degree_next_round;
	vector<int> right_degree_next_round;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
			left_vertices_to_be_peeled.push_back(u);
		}
	}
	int right_remain_nodes_num = g.num_v2;
	vector<vid_t> right_remain_nodes; right_remain_nodes.resize(g.num_v2);
	for (int i = 0; i < right_remain_nodes.size(); i++) {
		right_remain_nodes[i] = i;
	}
	int right_remain_nodes_tmp_num = 0;
	vector<vid_t> right_remain_nodes_tmp; right_remain_nodes_tmp.resize(g.num_v2);
	bool update_flag = false;
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		if (right_k - 1 > 0) {
			update_flag = true;
		}
		int pre_ = right_k - 1;
		bool stop = true;
		right_remain_nodes_tmp_num = 0;
		for (int i = 0; i < right_remain_nodes_num; i++) {
			vid_t v = right_remain_nodes[i];
			if (!g.right_delete[v]) {
				stop = false;
				right_remain_nodes_tmp[right_remain_nodes_tmp_num] = v;
				right_remain_nodes_tmp_num++;
				if (g.degree_v2[v] < right_k) {
					right_vertices_to_be_peeled.push_back(v);
				}
			}
		}
		swap(right_remain_nodes, right_remain_nodes_tmp);
		right_remain_nodes_num = right_remain_nodes_tmp_num;
		if (stop) break;
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			int oo_ = left_vertices_to_be_peeled.size();
			for (int j = 0; j < oo_; j++) {
				vid_t u = left_vertices_to_be_peeled[j];
				if (g.left_delete[u]) continue;
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int k = 0; k < ss; k++) {
					vid_t v = tmp_neigh_[k];
					if (g.right_delete[v]) continue;
					dd_ = --g.degree_v2[v];
					if (update_flag && dd_ == 0) {
						crossUpdate_for_kcore(g, left_k, pre_, v);
						g.right_delete[v] = true;
					}
					if (dd_ == pre_) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
				g.degree_v1[u] = 0;
				g.left_delete[u] = true;
				if (update_flag) {
					g.left_index[u][left_k] = pre_;
				}
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			oo_ = right_vertices_to_be_peeled.size();
			for (int j = 0; j < oo_; j++) {
				vid_t v = right_vertices_to_be_peeled[j];
				if (g.right_delete[v]) continue;
				vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v];
				int ss = tmp_neigh_.size();
				for (int k = 0; k < ss; k++) {
					vid_t u = tmp_neigh_[k];
					if (g.left_delete[u]) continue;
					dd_ = --g.degree_v1[u];
					if (update_flag && dd_ == 0) {
						g.left_index[u][left_k] = pre_;
						g.left_delete[u] = true;
					}
					if (dd_ == pre_left_k_) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.degree_v2[v] = 0;
				g.right_delete[v] = true;
				if (update_flag) {
					crossUpdate_for_kcore(g, left_k, pre_, v);
				}
			}
			right_vertices_to_be_peeled.clear();
		}
		if (right_k == 1) {
			left_degree_next_round = g.degree_v1;
			right_degree_next_round = g.degree_v2;
			left_deletion_next_round = g.left_delete;
			right_deletion_next_round = g.right_delete;
		}
	}
	//cout << 222 << endl;
	g.degree_v1 = left_degree_next_round;
	g.degree_v2 = right_degree_next_round;
	g.left_delete = left_deletion_next_round;
	g.right_delete = right_deletion_next_round;
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}

void alphaCopyPeel_for_kcore_with_pruning_rules(int left_k, BiGraph& g, bool skip, int& skip_num) {
	if (skip) {
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (!g.left_delete[u]) {
				g.left_index[u][left_k] = g.left_index[u][left_k - 1];
			}
		}
		for (vid_t v = 0; v < g.num_v2; v++) {
			if (!g.right_delete[v]) {
				int beta = -1;
				for (beta = g.right_index[v].size() - 1; beta > 0; beta--) {
					if (g.right_index[v][beta] == left_k - 1) {
						break;
					}
				}
				for (int x = beta; x > 0; x--) {
					if (g.right_index[v][x] < left_k) {
						g.right_index[v][x] = left_k;
					}
					else {
						break;
					}
				}
			}
		}
		return;
	}
	int min_left_degree_all_state = g.v1_max_degree;
	int dd_;
	int pre_left_k_ = left_k - 1;
	vector<bool> left_deletion_next_round;
	vector<bool> right_deletion_next_round;
	vector<int> left_degree_next_round;
	vector<int> right_degree_next_round;
	vector<vid_t> left_vertices_to_be_peeled;
	vector<vid_t> right_vertices_to_be_peeled;
	for (vid_t u = 0; u < g.getV1Num(); u++) {
		if (g.degree_v1[u] < left_k && !g.left_delete[u]) {
			left_vertices_to_be_peeled.push_back(u);
		}
	}
	int right_remain_nodes_num = g.num_v2;
	vector<vid_t> right_remain_nodes; right_remain_nodes.resize(g.num_v2);
	for (int i = 0; i < right_remain_nodes.size(); i++) {
		right_remain_nodes[i] = i;
	}
	int right_remain_nodes_tmp_num = 0;
	vector<vid_t> right_remain_nodes_tmp; right_remain_nodes_tmp.resize(g.num_v2);
	bool update_flag = false;
	for (int right_k = 1; right_k <= g.v2_max_degree + 1; right_k++) {
		if (right_k - 1 > 0) {
			update_flag = true;
		}
		int pre_ = right_k - 1;
		bool stop = true;
		right_remain_nodes_tmp_num = 0;
		for (int i = 0; i < right_remain_nodes_num; i++) {
			vid_t v = right_remain_nodes[i];
			if (!g.right_delete[v]) {
				stop = false;
				right_remain_nodes_tmp[right_remain_nodes_tmp_num] = v;
				right_remain_nodes_tmp_num++;
				if (g.degree_v2[v] < right_k) {
					right_vertices_to_be_peeled.push_back(v);
				}
			}
		}
		swap(right_remain_nodes, right_remain_nodes_tmp);
		right_remain_nodes_num = right_remain_nodes_tmp_num;
		if (stop) break;
		while (!left_vertices_to_be_peeled.empty() || !right_vertices_to_be_peeled.empty()) {
			// peel left
			int oo_ = left_vertices_to_be_peeled.size();
			for (int j = 0; j < oo_; j++) {
				vid_t u = left_vertices_to_be_peeled[j];
				if (g.left_delete[u]) continue;
				vector<vid_t>& tmp_neigh_ = g.neighbor_v1[u];
				int ss = tmp_neigh_.size();
				for (int k = 0; k < ss; k++) {
					vid_t v = tmp_neigh_[k];
					if (g.right_delete[v]) continue;
					dd_ = --g.degree_v2[v];
					if (update_flag && dd_ == 0) {
						crossUpdate_for_kcore(g, left_k, pre_, v);
						g.right_delete[v] = true;
					}
					if (dd_ == pre_) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
				g.degree_v1[u] = 0;
				g.left_delete[u] = true;
				if (update_flag) {
					g.left_index[u][left_k] = pre_;
				}
			}
			left_vertices_to_be_peeled.clear();
			// peel right
			oo_ = right_vertices_to_be_peeled.size();
			for (int j = 0; j < oo_; j++) {
				vid_t v = right_vertices_to_be_peeled[j];
				if (g.right_delete[v]) continue;
				vector<vid_t>& tmp_neigh_ = g.neighbor_v2[v];
				int ss = tmp_neigh_.size();
				for (int k = 0; k < ss; k++) {
					vid_t u = tmp_neigh_[k];
					if (g.left_delete[u]) continue;
					dd_ = --g.degree_v1[u];
					if (update_flag && dd_ == 0) {
						g.left_index[u][left_k] = pre_;
						g.left_delete[u] = true;
					}
					if (dd_ == pre_left_k_) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.degree_v2[v] = 0;
				g.right_delete[v] = true;
				if (update_flag) {
					crossUpdate_for_kcore(g, left_k, pre_, v);
				}
			}
			right_vertices_to_be_peeled.clear();
		}
		if (right_k == 1) {
			left_degree_next_round = g.degree_v1;
			right_degree_next_round = g.degree_v2;
			left_deletion_next_round = g.left_delete;
			right_deletion_next_round = g.right_delete;
		}
		if (min_left_degree_all_state > left_k) {
			for (vid_t u = 0; u < g.num_v1; u++) {
				if (!g.left_delete[u]) {
					if (g.degree_v1[u] < min_left_degree_all_state) {
						min_left_degree_all_state = g.degree_v1[u];
					}
				}
			}
		}
	}
	skip_num = min_left_degree_all_state;
	g.degree_v1 = left_degree_next_round;
	g.degree_v2 = right_degree_next_round;
	g.left_delete = left_deletion_next_round;
	g.right_delete = right_deletion_next_round;
	g.v1_max_degree = 0;
	g.v2_max_degree = 0;
	for (vid_t u = 0; u < g.degree_v1.size(); u++) {
		if (g.v1_max_degree < g.degree_v1[u]) g.v1_max_degree = g.degree_v1[u];
	}
	for (vid_t v = 0; v < g.degree_v2.size(); v++) {
		if (g.v2_max_degree < g.degree_v2[v]) g.v2_max_degree = g.degree_v2[v];
	}
}

void coreIndexBasic(BiGraph& g) {
	int left_degree_max = 0;
	int m = 0;
	for (int i = 0; i < g.getV1Num(); i++) {
		if (left_degree_max < g.getV1Degree(i)) left_degree_max = g.getV1Degree(i);
		m += g.getV1Degree(i);
	}
	int right_degree_max = 0;
	for (int i = 0; i < g.getV2Num(); i++) {
		if (right_degree_max < g.getV2Degree(i)) right_degree_max = g.getV2Degree(i);
	}
	// init g's max degree and index
	g.v1_max_degree = left_degree_max;
	g.v2_max_degree = right_degree_max;
	g.left_index.resize(g.getV1Num());
	g.right_index.resize(g.getV2Num());
	g.left_delete.resize(g.getV1Num());
	g.right_delete.resize(g.getV2Num());
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
		__alphaCopyPeel(left_k, g);
		if (PRINT_INDEX_DETAIL) {
			print_index_detail(g);
			cout << endl;
		}
	}
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = left_degree_max;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = right_degree_max;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}

	if (PRINT_INDEX_DETAIL) {
		cout << "inverse" << endl;
		cout << endl;
	}
	inv(g);
	for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
		__alphaCopyPeel(left_k, g);
		if (PRINT_INDEX_DETAIL) {
			inv(g);
			print_index_detail(g);
			cout << endl;
			inv(g);
		}
	}
	inv(g);
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = left_degree_max;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = right_degree_max;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}
}

void __coreIndexBasic(BiGraph& g) {
	int left_degree_max = 0;
	for (int i = 0; i < g.getV1Num(); i++) {
		if (left_degree_max < g.getV1Degree(i)) left_degree_max = g.getV1Degree(i);
	}
	int right_degree_max = 0;
	for (int i = 0; i < g.getV2Num(); i++) {
		if (right_degree_max < g.getV2Degree(i)) right_degree_max = g.getV2Degree(i);
	}
	// init g's max degree and index
	g.v1_max_degree = left_degree_max;
	g.v2_max_degree = right_degree_max;
	g.left_index.resize(g.getV1Num());
	g.right_index.resize(g.getV2Num());
	g.left_delete.resize(g.getV1Num());
	g.right_delete.resize(g.getV2Num());
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	BiGraph tmp;
	for (int i = 1; i <= g.v1_max_degree; i++) {
		tmp = g;
		alphaPeel(i, tmp, g);
	}
	// calculate right
	inv(g);
	for (int i = 1; i <= g.v1_max_degree; i++) {
		tmp = g;
		alphaPeel(i, tmp, g);
	}
	inv(g);
}

int coreIndexKCore(BiGraph& g) {
	int left_degree_max = 0;
	for (int i = 0; i < g.getV1Num(); i++) {
		if (left_degree_max < g.getV1Degree(i)) left_degree_max = g.getV1Degree(i);
	}
	int right_degree_max = 0;
	for (int i = 0; i < g.getV2Num(); i++) {
		if (right_degree_max < g.getV2Degree(i)) right_degree_max = g.getV2Degree(i);
	}
	// init g's max degree and index
	g.v1_max_degree = left_degree_max;
	g.v2_max_degree = right_degree_max;
	g.left_index.resize(g.getV1Num());
	g.right_index.resize(g.getV2Num());
	g.left_delete.resize(g.getV1Num());
	g.right_delete.resize(g.getV2Num());
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);//参数包括 : 一个迭代器，一个计数器以及一个值。该函数从迭代器指向的元素开始，将指定数量的元素设置为给定的值。
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	//int km = maxkcorenumber_test_optimal(g);
	//cout << "max kcore number: " << km << endl;
	int beta_s = 0;
	for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
		alphaCopyPeel_for_kcore(left_k, g);
		if (PRINT_INDEX_DETAIL) {
			print_index_detail(g);
			cout << endl;
		}
		beta_s = 0;
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.degree_v1[u] <= left_k) continue;
			int right_k = g.left_index[u][left_k];
			if (beta_s < right_k) beta_s = right_k;
		}
		if (beta_s <= left_k) break;
	}
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = left_degree_max;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = right_degree_max;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}

	if (PRINT_INDEX_DETAIL) {
		cout << "inverse" << endl;
		cout << endl;
	}
	inv(g);
	for (int left_k = 1; left_k <= beta_s; left_k++) {
		alphaCopyPeel_for_kcore(left_k, g);
		if (PRINT_INDEX_DETAIL) {
			inv(g);
			print_index_detail(g);
			cout << endl;
			inv(g);
		}
	}
	inv(g);
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = left_degree_max;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = right_degree_max;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}
	return beta_s;
}


int coreIndexKCore_with_pruning_rules(BiGraph& g) {
	int skip_times = 0;
	int left_degree_max = 0;
	for (int i = 0; i < g.getV1Num(); i++) {
		if (left_degree_max < g.getV1Degree(i)) left_degree_max = g.getV1Degree(i);
	}
	int right_degree_max = 0;
	for (int i = 0; i < g.getV2Num(); i++) {
		if (right_degree_max < g.getV2Degree(i)) right_degree_max = g.getV2Degree(i);
	}
	// init g's max degree and index
	g.v1_max_degree = left_degree_max;
	g.v2_max_degree = right_degree_max;
	g.left_index.resize(g.getV1Num());
	g.right_index.resize(g.getV2Num());
	g.left_delete.resize(g.getV1Num());
	g.right_delete.resize(g.getV2Num());
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	for (int i = 0; i < g.getV1Num(); i++) {
		g.left_index[i].resize(g.getV1Degree(i) + 1);
		fill_n(g.left_index[i].begin(), g.left_index[i].size(), 0);
	}
	for (int i = 0; i < g.getV2Num(); i++) {
		g.right_index[i].resize(g.getV2Degree(i) + 1);
		fill_n(g.right_index[i].begin(), g.right_index[i].size(), 0);
	}
	//int km = maxkcorenumber_test_optimal(g);
	//cout << "max kcore number: " << km << endl;
	int beta_s = 0;
	int skip_num = -1;
	for (int left_k = 1; left_k <= g.v1_max_degree; left_k++) {
		//cout << skip_num << endl;
		if (skip_num < left_k) {
			alphaCopyPeel_for_kcore_with_pruning_rules(left_k, g, false, skip_num);
		}
		else {
			alphaCopyPeel_for_kcore_with_pruning_rules(left_k, g, true, skip_num);
			skip_times++;
		}
		if (PRINT_INDEX_DETAIL) {
			print_index_detail(g);
			cout << endl;
		}
		beta_s = 0;
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (g.degree_v1[u] <= left_k) continue;
			int right_k = g.left_index[u][left_k];
			if (beta_s < right_k) beta_s = right_k;
		}
		if (beta_s <= left_k) break;
	}
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = left_degree_max;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = right_degree_max;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}

	if (PRINT_INDEX_DETAIL) {
		cout << "inverse" << endl;
		cout << endl;
	}
	inv(g);
	skip_num = -1;
	for (int left_k = 1; left_k <= beta_s; left_k++) {
		if (skip_num < left_k) {
			alphaCopyPeel_for_kcore_with_pruning_rules(left_k, g, false, skip_num);
		}
		else {
			alphaCopyPeel_for_kcore_with_pruning_rules(left_k, g, true, skip_num);
			skip_times++;
		}
		if (PRINT_INDEX_DETAIL) {
			inv(g);
			print_index_detail(g);
			cout << endl;
			inv(g);
		}
	}
	inv(g);
	// restore g
	fill_n(g.left_delete.begin(), g.left_delete.size(), false);
	g.v1_max_degree = left_degree_max;
	for (vid_t u = 0; u < g.num_v1; u++) {
		g.degree_v1[u] = g.neighbor_v1[u].size();
	}
	fill_n(g.right_delete.begin(), g.right_delete.size(), false);
	g.v2_max_degree = right_degree_max;
	for (vid_t v = 0; v < g.num_v2; v++) {
		g.degree_v2[v] = g.neighbor_v2[v].size();
	}
	cout << "skip_times: " << skip_times << endl;
	return beta_s;
}

void findabcore_core_index(int alpha, int beta, BiGraph& g, vector<bool>& left_sub, vector<bool>& right_sub) {
	for (vid_t u = 0; u < g.num_v1; u++) {
		if (g.degree_v1[u] < alpha) continue;
		if (g.left_index[u][alpha] >= beta) left_sub[u] = true;
	}
	for (vid_t v = 0; v < g.num_v2; v++) {
		if (g.degree_v2[v] < beta) continue;
		if (g.right_index[v][beta] >= alpha) right_sub[v] = true;
	}
}

bool test_abcore_equal(vector<bool>& left1, vector<bool>& right1, vector<bool>& left2, vector<bool>& right2) {
	for (int i = 0; i < left1.size(); i++) {
		if (left1[i] != left2[i]) {
			return false;
		}
	}
	for (int i = 0; i < right1.size(); i++) {
		if (right1[i] != right2[i]) {
			return false;
		}
	}
	return true;
}