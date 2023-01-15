#include "dyn_rebuild.h"
#include <cmath>
#include <stack>
using namespace std;

namespace Dyn_rebuild
{
	mutex g_mutex;
	volatile int alpha_mutex(0);
	volatile int beta_mutex(0);


	int calculate_ualpha_with_fixed_beta(BiGraph& g, vid_t u, int beta) {
		int ss = g.left_index[u].size();
		for (int alpha = 1; alpha < ss; alpha++) {
			if (g.left_index[u][alpha] < beta) {
				return alpha - 1;
			}
		}
		return ss - 1;
	}
	int calculate_vbeta_with_fixed_alpha(BiGraph& g, vid_t v, int alpha) {
		int ss = g.right_index[v].size();
		for (int beta = 1; beta < ss; beta++) {
			if (g.right_index[v][beta] < alpha) {
				return beta - 1;
			}
		}
		return ss - 1;
	}

	void compute_a_b_core(BiGraph& g, int alpha, int beta) {
		// intialize g.left_delete and g.right_deletion
		//fill_n(g.left_delete.begin(), g.left_delete.size(), false);
		//fill_n(g.right_delete.begin(), g.right_delete.size(), false);
		vector<vid_t> left_vertices_to_be_peeled;
		vector<vid_t> right_vertices_to_be_peeled;
		bool stop = true;
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (!g.left_delete[u]) {
				if (g.degree_v1[u] < alpha) {
					left_vertices_to_be_peeled.push_back(u);
				}
			}
		}
		for (vid_t v = 0; v < g.num_v2; v++) {
			if (!g.right_delete[v]) {
				if (g.degree_v2[v] < beta) {
					right_vertices_to_be_peeled.push_back(v);
				}
			}
		}
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
					int dd_ = --g.degree_v2[v];
					if (dd_ == 0) {
						// core part
						g.right_delete[v] = true;
					}
					if (dd_ == beta - 1) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
				g.degree_v1[u] = 0;
				g.left_delete[u] = true;
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
					int dd_ = --g.degree_v1[u];
					if (dd_ == 0) {
						g.left_delete[u] = true;
					}
					if (dd_ == alpha - 1) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.degree_v2[v] = 0;
				g.right_delete[v] = true;
			}
			right_vertices_to_be_peeled.clear();
		}
	}
	void compute_a_b_core_unchange(BiGraph& g, int alpha, int beta) {
		// intialize g.left_delete and g.right_deletion
		//fill_n(g.left_delete.begin(), g.left_delete.size(), false);
		//fill_n(g.right_delete.begin(), g.right_delete.size(), false);
		vector<vid_t> left_vertices_to_be_peeled;
		vector<vid_t> right_vertices_to_be_peeled;
		bool stop = true;
		for (vid_t u = 0; u < g.num_v1; u++) {
			if (!g.left_delete[u]) {
				if (g.degree_v1[u] < alpha) {
					left_vertices_to_be_peeled.push_back(u);
				}
			}
		}
		for (vid_t v = 0; v < g.num_v2; v++) {
			if (!g.right_delete[v]) {
				if (g.degree_v2[v] < beta) {
					right_vertices_to_be_peeled.push_back(v);
				}
			}
		}
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
					int dd_ = --g.degree_v2[v];
					if (dd_ == 0) {
						// core part
						g.right_delete[v] = true;
					}
					if (dd_ == beta - 1) {
						right_vertices_to_be_peeled.push_back(v);
					}
				}
				g.degree_v1[u] = 0;
				g.left_delete[u] = true;
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
					int dd_ = --g.degree_v1[u];
					if (dd_ == 0) {
						g.left_delete[u] = true;
					}
					if (dd_ == alpha - 1) {
						left_vertices_to_be_peeled.push_back(u);
					}
				}
				g.degree_v2[v] = 0;
				g.right_delete[v] = true;
			}
			right_vertices_to_be_peeled.clear();
		}
		g.degree_v1 = g.degree_of_left_node_initial;
		g.degree_v2 = g.degree_of_right_node_initial;
	}

	int upgrade_vbeta_insertion(BiGraph& g, vid_t v, int alpha, int beta) {
		int a1 = g.right_index[v][beta];
		if (beta == 0) {
			a1 = alpha;
		}
		int a2 = g.right_index[v][beta + 1];
		if (a1 < alpha) {
			return -1;
		}
		else if (a1 >= alpha && a2 < alpha) {
			return beta;
		}
		else {
			return 1e9;
		}
	}

	int upgrade_vbeta_removal(BiGraph& g, vid_t v, int alpha, int beta) {
		int a1 = g.right_index[v][beta];
		int a2 = -1;
		if (beta + 1 < g.right_index[v].size()) {
			a2 = g.right_index[v][beta + 1];
		}
		if (a1 < alpha) {
			return -1;
		}
		else if (a1 >= alpha && a2 < alpha) {
			return beta;
		}
		else {
			return 1e9;
		}
	}
	// may incur errors !!!
	int compute_vbeta_regarding_specific_alpha(BiGraph& g, vid_t v, int alpha) {
		int vbeta = -1;
		bool changed = false;
		for (int b = 1; b <= g.right_index[v].size() - 1; b++) {
			if (g.right_index[v][b] < alpha) {
				vbeta = b - 1;
				changed = true;
				break;
			}
		}
		if (!changed) {
			vbeta = g.right_index[v].size() - 1;
		}
		if (vbeta < 0) {
			//cout << "error: vbeta" << endl;
		}
		return vbeta;
	}
	int compute_vbeta_regarding_specific_alpha__(BiGraph& g, vid_t v, int alpha) {
		int vbeta = -1;
		bool changed = false;
		for (int b = 1; b <= g.degree_v2[v]; b++) {
			if (g.right_index[v][b] < alpha) {
				vbeta = b - 1;
				changed = true;
				break;
			}
		}
		if (!changed) {
			vbeta = g.degree_v2[v];
		}
		if (vbeta < 0) {
			cout << "error: vbeta" << endl;
		}
		return vbeta;
	}

	void casacade_remove_candidate_nodes(vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right) {
		if (left_side) {
			candidate_nodes_left.erase(node);
			for (auto v : visited_nodes_2_info_left[node].neighbors_in_candidates) {
				if (candidate_nodes_right.find(v) != candidate_nodes_right.end()) {
					visited_nodes_2_info_right[v].current_supports--;
					if (visited_nodes_2_info_right[v].current_supports <= beta) {
						casacade_remove_candidate_nodes(v, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
		else {
			candidate_nodes_right.erase(node);
			for (auto u : visited_nodes_2_info_right[node].neighbors_in_candidates) {
				if (candidate_nodes_left.find(u) != candidate_nodes_left.end()) {
					visited_nodes_2_info_left[u].current_supports--;
					if (visited_nodes_2_info_left[u].current_supports < alpha) {
						casacade_remove_candidate_nodes(u, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
	}

	void casacade_insert_candidate_nodes(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_right) {
		if (left_side) {
			candidate_nodes_left.insert(node);
			for (auto v : visited_nodes_2_info_left[node].neighbors_not_in_candidates) {
				if (candidate_nodes_right.find(v) == candidate_nodes_right.end()) {
					visited_nodes_2_info_right[v].current_supports--;
					if (visited_nodes_2_info_right[v].current_supports < beta) {
						casacade_insert_candidate_nodes(g, v, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
			for (auto v : visited_nodes_2_info_left[node].nodes_to_be_visited) {
				if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
					dyn_update_removal_dfs(g, v, alpha, beta, !left_side,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
		}
		else {
			candidate_nodes_right.insert(node);
			for (auto u : visited_nodes_2_info_right[node].neighbors_not_in_candidates) {
				if (candidate_nodes_left.find(u) == candidate_nodes_left.end()) {
					visited_nodes_2_info_left[u].current_supports--;
					if (visited_nodes_2_info_left[u].current_supports < alpha) {
						casacade_insert_candidate_nodes(g, u, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
			for (auto u : visited_nodes_2_info_right[node].nodes_to_be_visited) {
				if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
					dyn_update_removal_dfs(g, u, alpha, beta, !left_side,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
		}
	}

	void dyn_update_removal_dfs(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_right) {
		if (!left_side && node == 4) {
			int kkkkkk = 0;// for debug
		}
		if (left_side) {
			visited_nodes_left.insert(node);
			vector<vid_t>& neigh = g.neighbor_v1[node];
			for (auto v : neigh) {
				// could be a candidate only if v's deg >= beta
				// aim at finding whether v still belongs to alpha-beta-core
				if (g.degree_v2[v] >= beta) {
					// could be optimized???
					int vbeta = upgrade_vbeta_removal(g, v, alpha, beta);
					if (vbeta > beta) {
						visited_nodes_2_info_left[node].current_supports++;
					}
					else if (vbeta == beta) {
						if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
							visited_nodes_2_info_left[node].current_supports++;
							visited_nodes_2_info_left[node].nodes_to_be_visited.push_back(v);
						}
						else if (visited_nodes_right.find(v) != visited_nodes_right.end() && candidate_nodes_right.find(v) == candidate_nodes_right.end()) {
							visited_nodes_2_info_left[node].current_supports++;
							visited_nodes_2_info_left[node].neighbors_not_in_candidates.push_back(v);
							visited_nodes_2_info_right[v].neighbors_not_in_candidates.push_back(node);
						}
						else {
							//do nothing
						}
					}
					else {
						// do nothing
					}
				}
			}
			// expand current node
			if (visited_nodes_2_info_left[node].current_supports < alpha) {
				casacade_insert_candidate_nodes(g, node, alpha, beta, left_side,
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
					visited_nodes_2_info_left, visited_nodes_2_info_right);
				for (auto v : visited_nodes_2_info_left[node].nodes_to_be_visited) {
					if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
						dyn_update_removal_dfs(g, v, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
		else {
			visited_nodes_right.insert(node);
			vector<vid_t>& neigh = g.neighbor_v2[node];
			for (auto u : neigh) {
				if (g.degree_v1[u] >= alpha) {
					int ubeta = g.left_index[u][alpha];
					if (ubeta > beta) {
						visited_nodes_2_info_right[node].current_supports++;
					}
					else if (ubeta == beta) {
						if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
							visited_nodes_2_info_right[node].current_supports++;
							visited_nodes_2_info_right[node].nodes_to_be_visited.push_back(u);
						}
						// can be simplified
						else if (visited_nodes_left.find(u) != visited_nodes_left.end() && candidate_nodes_left.find(u) == candidate_nodes_left.end()) {
							visited_nodes_2_info_right[node].current_supports++;
							// this node is doomed! must be the end point of insertion edge otherwise it can not be visited!
							if (g.right_index[node][beta] < alpha || g.degree_v2[node] < beta) continue;
							visited_nodes_2_info_right[node].neighbors_not_in_candidates.push_back(u);
							visited_nodes_2_info_left[u].neighbors_not_in_candidates.push_back(node);
						}
						else {
							//do nothing
						}
					}
					else {
						// do nothing
					}
				}
			}
			if (visited_nodes_2_info_right[node].current_supports < beta) {
				casacade_insert_candidate_nodes(g, node, alpha, beta, left_side,
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
					visited_nodes_2_info_left, visited_nodes_2_info_right);
				for (auto u : visited_nodes_2_info_right[node].nodes_to_be_visited) {
					if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
						dyn_update_removal_dfs(g, u, alpha, beta, !left_side,
							visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
							visited_nodes_2_info_left, visited_nodes_2_info_right);
					}
				}
			}
		}
	}

	int compute_ualpha_regarding_specific_beta(BiGraph& g, vid_t u, int beta) {
		int ualpha = -1;
		bool changed = false;
		for (int b = 1; b <= g.left_index[u].size() - 1; b++) {
			if (g.left_index[u][b] < beta) {
				ualpha = b - 1;
				changed = true;
				break;
			}
		}
		if (!changed) {
			ualpha = g.left_index[u].size() - 1;
		}
		if (ualpha < 0) {
			//cout << "error: vbeta" << endl;
		}
		return ualpha;
	}

	int comput_u_alpha(std::vector<skyline_block*>& skyline_index_uv, std::vector<skyline_block*>& skyline_index_uv_reverse, vid_t uv, int alpha_or_beta, int num){
        if(alpha_or_beta == 0) return 0;
        if(skyline_index_uv_reverse[uv]->mapset.count(alpha_or_beta)!=0) return skyline_index_uv_reverse[uv]->mapset[alpha_or_beta];
        if(alpha_or_beta > skyline_index_uv[uv]->mapset.begin()->second) return 0;
        if(alpha_or_beta < skyline_index_uv_reverse[uv]->mapset.begin()->first) return skyline_index_uv_reverse[uv]->mapset.begin()->second;
        for (auto iter = skyline_index_uv[uv]->mapset.begin(); iter != skyline_index_uv[uv]->mapset.end(); iter++){
            if (alpha_or_beta > iter->second) return (--iter)->first;
        }
	}

	int find_kth(vector<int> B, int n, int k){
        if (k < 0||k > n) return 0;
        sort(B.begin(), B.end());
        return B[n - k];
	}

	int find_skyline(std::vector<skyline_block*>& skyline_index_uv, vid_t uv, int alpha_or_beta, int num){
	    auto start_find = chrono::system_clock::now();
        if (skyline_index_uv[uv]->mapset.count(alpha_or_beta) != 0) {
                return skyline_index_uv[uv]->mapset[alpha_or_beta];
        }
        else{
            if (alpha_or_beta <= (--skyline_index_uv[uv]->mapset.end())->first||alpha_or_beta <= skyline_index_uv[uv]->mapset.begin()->first){
                int i = alpha_or_beta;
                while (skyline_index_uv[uv]->mapset.count(i) == 0 && i <= num && (i <= (--skyline_index_uv[uv]->mapset.end())->first || i <= skyline_index_uv[uv]->mapset.begin()->first)) {
                    i++;
                }
                return skyline_index_uv[uv]->mapset[i];
            }
            else {
                return 0;
            }
        }
	}

	void build_skyline(BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_u_reverse, 
	std::vector<skyline_block*>& skyline_index_v, std::vector<skyline_block*>& skyline_index_v_reverse){
		for (int i = 0; i < g.left_index_delete.size(); i++) {
			for (int j = 1; j < g.left_index_delete[i].size(); j++) {
				if(g.left_index_delete[i][j] == 0) {
					if(g.left_index_change[i][j] !=0){
						skyline_index_u[i]->mapset[j] = g.left_index_change[i][j];
						skyline_index_u_reverse[i]->mapset[g.left_index_change[i][j]] = j;
						g.left_index[i][j] = g.left_index_change[i][j];
					}
					else{
						skyline_index_u[i]->mapset.erase(j);
						skyline_index_u_reverse[i]->mapset.erase(g.left_index_change[i][j]);
						g.left_index[i].pop_back();
					}
				}
			}
		}

		for (int i = 0; i < g.right_index_delete.size(); i++) {
			for (int j = 1; j < g.right_index_delete[i].size(); j++) {
				if(g.right_index_delete[i][j] == 0){
					if(g.right_index_change[i][j] != 0){
						skyline_index_v[i] -> mapset[g.right_index_change[i][j]] = j;
						skyline_index_v_reverse[i] -> mapset[j] = g.right_index_change[i][j];
						g.right_index[i][j] = g.right_index_change[i][j];
					}
					else{
						skyline_index_v[i] -> mapset.erase(g.right_index_change[i][j]);
						skyline_index_v_reverse[i] -> mapset.erase(j);
						g.right_index[i].pop_back();
					}
				}
			}
		}

		for (int i = 0; i < g.getV1Num(); i++) {
			fill_n(g.left_index_change[i].begin(), g.left_index_change[i].size(), 0);
		}
		for (int i = 0; i < g.getV2Num(); i++) {
			fill_n(g.right_index_change[i].begin(), g.right_index_change[i].size(), 0);
		}

		g.left_index_delete = g.left_index;
		g.right_index_delete = g.right_index;
	}

	double update_skyline_index_swap_with_bfs_delete(BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v,
                                                  std::vector<skyline_block*>& skyline_index_u_reverse, std::vector<skyline_block*>& skyline_index_v_reverse, vid_t u, vid_t v){
		auto start = chrono::system_clock::now();
		g.deleteEdge(u, v);

        map<int,int> pair_to_be_visit;
        map<int,int> pair_to_be_visit_v;
        for(auto it : skyline_index_u[u]->mapset){
            int alpha = it.first; int beta = it.second;
            if(find_skyline(skyline_index_v, v, alpha, g.num_v2) >= beta){
			//if(compute_vbeta_regarding_specific_alpha(g, v, alpha) >= beta){
                pair_to_be_visit[alpha] = beta;
            }
        }
        for(auto it : skyline_index_v[v]->mapset){
            int alpha = it.first; int beta = it.second;
            if(find_skyline(skyline_index_u, u, alpha,g.num_v1) >= beta){
			//if(g.left_index[u][alpha] >= beta){
                //if(pair_to_be_visit.count(alpha) != 0){
                //    if (pair_to_be_visit[alpha] == beta) continue;
                //}
                pair_to_be_visit_v[alpha] = beta;
            }
        }

        if(pair_to_be_visit.size() == 0){
            pair_to_be_visit[g.degree_v1[u] + 1] = skyline_index_u[u] ->mapset[g.degree_v1[u] + 1];
        }
        else{
            if(pair_to_be_visit.count(g.degree_v1[u] + 1) == 0){
                pair_to_be_visit[g.degree_v1[u] + 1] = find_skyline(skyline_index_u, u, g.degree_v1[u] + 1,g.num_v1);
            }
        }
        if(pair_to_be_visit_v.size() == 0){
            pair_to_be_visit_v[skyline_index_u_reverse[u] ->mapset[g.degree_v2[v] + 1]] = g.degree_v2[v] + 1;
        }
        else{
            if(pair_to_be_visit_v.count(g.degree_v2[v] + 1) == 0){
                pair_to_be_visit_v[find_skyline(skyline_index_v_reverse, v, g.degree_v2[v] + 1,g.num_v2)] = g.degree_v2[v] + 1;
            }
        }

		unordered_set<vid_t> visited_nodes_left; unordered_set<vid_t> visited_nodes_right;
		unordered_set<vid_t> candidate_nodes_left; unordered_set<vid_t> candidate_nodes_right;
		unordered_map<vid_t, traversal_info_removal> visited_nodes_2_info_left; unordered_map<vid_t, traversal_info_removal> visited_nodes_2_info_right;
        
        vector<int> vector_A;
		vector<vid_t>& neigh_1 = g.neighbor_v1[u];
		vector<vid_t>& neigh_2 = g.neighbor_v2[v];
        for (auto iter = pair_to_be_visit.rbegin(); iter != pair_to_be_visit.rend(); iter++){
            visited_nodes_left.clear(); visited_nodes_right.clear();
            candidate_nodes_left.clear(); candidate_nodes_right.clear();
            visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
            int current_alpha = iter->first;
            int current_beta = iter->second;
			if (current_alpha == 0 || current_beta == 0) continue;
			dyn_update_removal_dfs(g, u, current_alpha, current_beta, true, visited_nodes_left, visited_nodes_right,
                                                                candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
			if(candidate_nodes_left.size() == 0 ) {
                pair_to_be_visit_v[current_alpha] = current_beta;
                continue;
            }
			if(pair_to_be_visit_v.count(current_alpha) != 0) {
                dyn_update_removal_dfs(g, v, current_alpha, current_beta, false, visited_nodes_left, visited_nodes_right,
                                                                candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
				pair_to_be_visit_v[current_alpha] = 0;
			}
            int beta_ = 0;
			int beta_2 = 0;
            vector<int> vector_B;
			vector<int> vector_C;
            if(current_alpha - 1 != 0){
                for (vid_t nbr_node:neigh_1){
                    //vector_B.push_back(find_skyline(skyline_index_v, nbr_node, current_alpha - 1,g.num_v2));
					vector_B.push_back(compute_vbeta_regarding_specific_alpha(g, nbr_node, current_alpha - 1));
					vector_C.push_back(compute_vbeta_regarding_specific_alpha(g, nbr_node, current_alpha));
                }
                beta_ = max(find_kth(vector_B,vector_B.size(),current_alpha - 1),current_beta);
                vector_B.clear();
				beta_2 = find_kth(vector_C,vector_C.size(),current_alpha);
            	vector_C.clear();
            }
			else{
				for (vid_t nbr_node:neigh_1){
					//vector_B.push_back(find_skyline(skyline_index_v, nbr_node, current_alpha,g.num_v2));
					vector_B.push_back(compute_vbeta_regarding_specific_alpha(g, nbr_node, current_alpha));
				}
				beta_2 = find_kth(vector_B,vector_B.size(),current_alpha);
				vector_B.clear();
			}

			for (auto w : candidate_nodes_left){
				g.left_index_delete[w][current_alpha] = 0;
				if (w == u){
					int i = current_alpha;
					while (g.left_index_change[w][i] < beta_2 && i < g.left_index[w].size()){
						if(beta_2 < g.left_index[w][i]) g.left_index_change[w][i] = beta_2;
						i++;
					}
				}
				else{
					int i = current_alpha;
					while (g.left_index_change[w][i] < current_beta - 1 && i < g.left_index[w].size()){
						if(current_beta - 1 < g.left_index[w][i]) {
							g.left_index_change[w][i] = current_beta - 1;
						}
						i++;
					}
				}
				if (g.left_index_change[w][current_alpha - 1] < beta_){
					if(beta_ < g.left_index[w][current_alpha - 1]) g.left_index_change[w][current_alpha - 1] = beta_;
				}
			}
			for (auto w : candidate_nodes_right) {
				g.right_index_delete[w][current_beta] = 0;
				if (g.right_index_change[w][current_beta - 1] < current_alpha){
					if(current_alpha < g.right_index[w][current_beta - 1]) g.right_index_change[w][current_beta - 1] = current_alpha;
				}
				if (w == v){
					int alpha_2 = 0;
					for (vid_t nbr_node:neigh_2){
						//vector_B.push_back(comput_u_alpha(skyline_index_u, skyline_index_u_reverse, nbr_node, current_beta, g.num_v1));
						vector_B.push_back(compute_ualpha_regarding_specific_beta(g, nbr_node, current_beta));
					}
					alpha_2 = find_kth(vector_B,vector_B.size(),current_beta);
					vector_B.clear();
					int i = current_beta;
					while (g.right_index_change[w][i] < alpha_2 && i < g.right_index[w].size()){
						if(alpha_2 < g.right_index[w][i]) g.right_index_change[w][i] = alpha_2;
						i++;
					}
				}
				else{
					for (int i = current_beta; i <= beta_; i++){ 
						if (g.right_index_change[w][i] < current_alpha - 1 && current_alpha - 1 < g.right_index[w][i]) g.right_index_change[w][i] = current_alpha - 1; 
					}
				}

			}
        }
        vector_A.clear();
        for (auto iter = pair_to_be_visit_v.begin(); iter != pair_to_be_visit_v.end(); iter++){
            int current_alpha = iter->first;
            int current_beta = iter->second;
			if (current_alpha == 0 || current_beta == 0) continue;
            //if(already_u.count(current_alpha) != 0) continue;
            visited_nodes_left.clear(); visited_nodes_right.clear();
            candidate_nodes_left.clear(); candidate_nodes_right.clear();
            visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
            dyn_update_removal_dfs(g, v, current_alpha, current_beta, false, visited_nodes_left, visited_nodes_right,
                                                                candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
			if(candidate_nodes_right.size() == 0 ) continue;
            int alpha_ = 0;
			int alpha_2 = 0;
			vector<int> vector_B;
			vector<int> vector_C;
            if(current_beta - 1 != 0){
                for (vid_t nbr_node:neigh_2){
                    //vector_B.push_back(comput_u_alpha(skyline_index_u, skyline_index_u_reverse, nbr_node, current_beta - 1, g.num_v1));
					vector_B.push_back(compute_ualpha_regarding_specific_beta(g, nbr_node, current_beta - 1));
					vector_C.push_back(compute_ualpha_regarding_specific_beta(g, nbr_node, current_beta));
                }
                alpha_ = max(find_kth(vector_B,vector_B.size(),current_beta - 1),current_alpha);
                vector_B.clear();
				alpha_2 = find_kth(vector_C,vector_C.size(),current_beta);
            	vector_C.clear();
            }
			else{
				for (vid_t nbr_node:neigh_2){
					//vector_B.push_back(comput_u_alpha(skyline_index_u, skyline_index_u_reverse, nbr_node, current_beta, g.num_v1));
					vector_B.push_back(compute_ualpha_regarding_specific_beta(g, nbr_node, current_beta));
				}
				alpha_2 = find_kth(vector_B,vector_B.size(),current_beta);
				vector_B.clear();
			}

			 
			for (auto w : candidate_nodes_right){
				g.right_index_delete[w][current_beta] = 0;
				if (w == v){
					int i = current_beta;
					while (g.right_index_change[w][i] < alpha_2 && i < g.right_index[w].size()){
						if(alpha_2 < g.right_index[w][i]) g.right_index_change[w][i] = alpha_2;
						i++;
					}
				}
				else{
					int i = current_beta;
					while (g.right_index_change[w][i] < current_alpha - 1 && i < g.right_index[w].size()){
						if(current_alpha - 1 < g.right_index[w][i]) g.right_index_change[w][i] = current_alpha - 1;
						i++;
					}
				}
				if (g.right_index_change[w][current_beta - 1] < alpha_){
					if(alpha_ < g.right_index[w][current_beta - 1])g.right_index_change[w][current_beta - 1] = alpha_;
				}
			}
			for (auto w : candidate_nodes_left) {
				g.left_index_delete[w][current_alpha] = 0;
				if (g.left_index_change[w][current_alpha - 1] < current_beta){
					if(current_beta < g.left_index[w][current_alpha - 1]) g.left_index_change[w][current_alpha - 1] = current_beta;
				}
				for (int i = current_alpha; i <= alpha_; i++){ 
					if (g.left_index_change[w][i] < current_beta - 1 && current_beta - 1 < g.left_index[w][i]) g.left_index_change[w][i] = current_beta - 1; 
				}
			}
        }
		g.left_index[u].pop_back();
		g.right_index[v].pop_back();

        auto end = chrono::system_clock::now();
		//build_skyline(g, skyline_index_u, skyline_index_u_reverse, skyline_index_v, skyline_index_v_reverse);
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
        //build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v, skyline_index_u_reverse, skyline_index_v_reverse);
        chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count();
    }
}
