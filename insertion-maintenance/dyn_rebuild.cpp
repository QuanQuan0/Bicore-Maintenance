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

	void dyn_update_insertion_dfs(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
		unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right,
		unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right) {
		// node + left_side * (g.num_v1 + g.num_v2) to distinguish between left and right
		vector<pair<pair<vid_t,vid_t>,bool>> stack; // pre cur left_side
		stack.push_back(make_pair(make_pair(-1,node), left_side));
		bool start_root = true;
		while (!stack.empty()) {
			auto back = stack.back();
			stack.pop_back();
			//vector<vid_t> nodes_to_be_visited;
			if (back.second) {
				vid_t pre = back.first.first;
				vid_t node = back.first.second;
				if (visited_nodes_left.find(node) != visited_nodes_left.end()) {
					continue;
				}
				if (candidate_nodes_right.find(pre) == candidate_nodes_right.end() && !start_root) {
					continue;
				}
				if (start_root) {
					start_root = false;
				}
				candidate_nodes_left.insert(node);
				visited_nodes_left.insert(node);
				vector<vid_t>& neigh = g.neighbor_v1[node];
				for (auto v : neigh) {
					// could be a candidate only if v's deg > beta
					// aim at finding whether v belongs to alpha-beta+1-core
					if (g.degree_v2[v] > beta) {
						// could be optimized???
						int vbeta = upgrade_vbeta_insertion(g, v, alpha, beta);
						//int vbeta = compute_vbeta_regarding_specific_alpha(g, v, alpha);
						if (vbeta > beta) {
							visited_nodes_2_info_left[node].current_supports++;
						}
						else if (vbeta == beta) {
							
							if (visited_nodes_right.find(v) == visited_nodes_right.end()) {
								visited_nodes_2_info_left[node].current_supports++;
								stack.push_back(make_pair(make_pair(node,v),false));
							}
							else if (visited_nodes_right.find(v) != visited_nodes_right.end() && candidate_nodes_right.find(v) != candidate_nodes_right.end()) {
								visited_nodes_2_info_left[node].current_supports++;
								visited_nodes_2_info_left[node].neighbors_in_candidates.push_back(v);
								visited_nodes_2_info_right[v].neighbors_in_candidates.push_back(node);
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
				if (visited_nodes_2_info_left[node].current_supports >= alpha) {
					// do nothing
				}
				else {
					casacade_remove_candidate_nodes(node, alpha, beta, true,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
			else {
				vid_t pre = back.first.first;
				vid_t node = back.first.second;
				if (visited_nodes_right.find(node) != visited_nodes_right.end()) {
					continue;
				}
				if (candidate_nodes_left.find(pre) == candidate_nodes_left.end() && !start_root) {
					continue;
				}
				if (start_root) {
					start_root = false;
				}
				candidate_nodes_right.insert(node);
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
								stack.push_back(make_pair(make_pair(node, u), true));
							}
							else if (visited_nodes_left.find(u) != visited_nodes_left.end() && candidate_nodes_left.find(u) != candidate_nodes_left.end()) {
								visited_nodes_2_info_right[node].current_supports++;
								visited_nodes_2_info_right[node].neighbors_in_candidates.push_back(u);
								visited_nodes_2_info_left[u].neighbors_in_candidates.push_back(node);
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
				if (visited_nodes_2_info_right[node].current_supports > beta) {
					
				}
				else {
					casacade_remove_candidate_nodes(node, alpha, beta, false,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}
			}
		}
	}

	int find_kth(vector<int> B, int n, int k){
		if (k < 0||k > n) return 0;
		sort(B.begin(), B.end());
		return B[n - k];
	}

	double update_skyline_index_swap_with_bfs(BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v, vid_t u, vid_t v){
		auto start = chrono::system_clock::now();
		g.addEdge(u, v);
		vector<vid_t>& nbr_u = g.neighbor_v1[u]; vector<vid_t>& nbr_v = g.neighbor_v2[v];
		if(g.degree_v1[u] == 1){
			g.left_index[u].push_back(g.degree_v2[v]);
			g.right_index[v].push_back(g.degree_v1[u]);
			for(auto nbr_node : nbr_v){
				if(g.left_index[nbr_node][1] < g.degree_v2[v]) g.left_index[nbr_node][1] = g.degree_v2[v];
			}
			auto end = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = end - start;
        	return elapsed_seconds.count();
		}
		if(g.degree_v2[v] == 1){
			g.right_index[v].push_back(g.degree_v1[u]);
			g.left_index[u].push_back(g.degree_v2[v]);
			for (auto nbr_node : nbr_u){
				if(g.right_index[nbr_node][1] < g.degree_v1[u]) g.right_index[nbr_node][1] = g.degree_v1[u];
			}
			auto end = chrono::system_clock::now();
			chrono::duration<double> elapsed_seconds = end - start;
        	return elapsed_seconds.count();
		}
		g.left_index[u].push_back(0);
		g.right_index[v].push_back(0);
		map<vector<int>, int> have_update;
        vector<int> B; vector<int> C;
		vector<int> D; vector<int> E;
		vector<int> F;
        vector<int> vector_A;
		unordered_set<vid_t> visited_nodes_left; unordered_set<vid_t> visited_nodes_right;
		unordered_set<vid_t> candidate_nodes_left; unordered_set<vid_t> candidate_nodes_right;
		unordered_map<vid_t, traversal_info_insertion> visited_nodes_2_info_left; unordered_map<vid_t, traversal_info_insertion> visited_nodes_2_info_right;
        
		for (auto nbr_node : nbr_u){
			B.push_back(compute_vbeta_regarding_specific_alpha(g, nbr_node, g.degree_v1[u]));
		}
		int beta_2 = find_kth(B, B.size(), g.degree_v1[u]);
		B.clear();
		visited_nodes_left.clear(); visited_nodes_right.clear();
		candidate_nodes_left.clear(); candidate_nodes_right.clear();
		visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
		g.left_index[u][g.degree_v1[u]] = max(beta_2,1);
		dyn_update_insertion_dfs(g, u, g.degree_v1[u], beta_2, true, visited_nodes_left, visited_nodes_right,
			candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		for (auto i : candidate_nodes_left){
			if(g.left_index[i][g.degree_v1[u]] < beta_2 + 1) {
				g.left_index[i][g.degree_v1[u]] = beta_2 + 1;
			}
		}
		for (auto i : candidate_nodes_right){
			if(g.right_index[i][beta_2 + 1] < g.degree_v1[u]) g.right_index[i][beta_2 + 1] = g.degree_v1[u];
		}

        for(int i_iter = skyline_index_u[u]->vector_set.size()-2; i_iter >= 0; i_iter-=2){
            int current_alpha = skyline_index_u[u]->vector_set[i_iter];
            int current_beta = skyline_index_u[u]->vector_set[i_iter+1];
			if(g.right_index[v].size() < current_beta) break;
			if (current_alpha > g.right_index[v][current_beta]){
				if (g.left_index[u][current_alpha+1] <= compute_vbeta_regarding_specific_alpha(g,v,current_alpha+1) && skyline_index_u[u]->vector_set[i_iter + 2] != current_alpha + 1){
					for (auto nbr_node : nbr_u){
						B.push_back(compute_vbeta_regarding_specific_alpha(g,nbr_node,current_alpha+1));
					}
					int beta = find_kth(B, B.size(),current_alpha + 1);
					if(beta >= skyline_index_u[u]->vector_set[i_iter+3]){
						if (g.left_index[u][current_alpha+1] < beta) {
							g.left_index[u][current_alpha+1] = beta;
							visited_nodes_left.clear(); visited_nodes_right.clear();
							candidate_nodes_left.clear(); candidate_nodes_right.clear();
							visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
							dyn_update_insertion_dfs(g, u, current_alpha+1, beta, true, visited_nodes_left, visited_nodes_right,
								candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
							for (auto i : candidate_nodes_left){
								g.left_index[i][current_alpha+1] = beta+1;
							}
							for (auto i : candidate_nodes_right){
								g.right_index[i][beta+1] = current_alpha+1;
							}
						}
					}
					B.clear();
				}
				continue;
			}
            if (current_beta > g.degree_v2[v]) continue;
			int value;
            for (auto nbr_node : nbr_u){
				value = compute_vbeta_regarding_specific_alpha(g, nbr_node, current_alpha + 1);
                B.push_back(compute_vbeta_regarding_specific_alpha(g, nbr_node, current_alpha));
				D.push_back(value);
				if (nbr_node != v) E.push_back(value);
            }
			int alpha_u_whether_alpha_add_ = find_kth(E, E.size(),current_alpha);
            E.clear();
			int beta_3 = find_kth(D, D.size(),current_alpha + 1);
			D.clear();
            int beta_ = find_kth(B, B.size(),current_alpha);
            if (beta_ == 0) beta_+= 1;
            int beta_v = compute_vbeta_regarding_specific_alpha(g, v, current_alpha);
            int beta_first = find_kth(B, B.size(),B.size() - current_alpha + 1);
            int beta_sub;
            if (current_alpha != 1) beta_sub = find_kth(B, B.size(),current_alpha - 1);
            else beta_sub = beta_v + 1;
            B.clear();
            for (auto nbr_node : nbr_v){
                B.push_back(compute_ualpha_regarding_specific_beta(g, nbr_node, beta_ + 1));
				if (nbr_node != u) D.push_back(g.left_index[nbr_node][current_alpha + 1]);
            }
            int alpha_whether_beta_add = find_kth(B, B.size(),beta_);
            B.clear();
            int alpha_u_whether_alpha_add;
            if (D.size() !=0 || current_beta - 1 != 0) alpha_u_whether_alpha_add = find_kth(D, D.size(),current_beta - 1);
            else alpha_u_whether_alpha_add = current_beta;
            D.clear();

            // u-alpha change
			if(current_alpha != g.degree_v1[u] -1){
				F.push_back(current_alpha + 1);
				F.push_back(beta_3 + 1);
				if(beta_v >= beta_3 && compute_vbeta_regarding_specific_alpha(g, v, current_alpha + 1) >= beta_3 - 1 && have_update.count(F) == 0){  
					if (alpha_u_whether_alpha_add >= beta_3 + 1 && alpha_u_whether_alpha_add_ >= beta_3 + 1){ 
						//skyline_index_u[u]->vector_set[i_iter] = current_alpha + 1;
						g.left_index[u][current_alpha + 1] = beta_3 + 1;
						C.push_back(current_alpha + 1);
						C.push_back(beta_3 + 1);
						g.is_node_change_in_pair_v1[u][C] = 1;
						C.clear();
						//have_update_u[current_alpha + 1] = current_beta;
					}
					else{
						g.left_index[u][current_alpha + 1] = beta_3;
						visited_nodes_left.clear(); visited_nodes_right.clear();
						candidate_nodes_left.clear(); candidate_nodes_right.clear();
						visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
						dyn_update_insertion_dfs(g, u, current_alpha + 1, beta_3, true, visited_nodes_left, visited_nodes_right,
							candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
						for (auto i : candidate_nodes_left){
							if(g.left_index[i][current_alpha + 1] < beta_3 + 1) g.left_index[i][current_alpha + 1] = beta_3 + 1;
						}
						for (auto i : candidate_nodes_right){
							if(g.right_index[i][beta_3 + 1] < current_alpha + 1) g.right_index[i][beta_3 + 1] = current_alpha + 1;
						}
					}
				}
				B.clear();
				if(g.left_index[u][current_alpha + 1] < beta_3) g.left_index[u][current_alpha + 1] = beta_3;
				if(g.left_index[u][current_alpha + 1] == beta_3 + 1){
					have_update[F] = 1;
					F.clear();
				}
			}
            
            // u-beta change
			F.push_back(current_alpha);
			F.push_back(beta_ + 1);
			if(have_update.count(F) != 0) {
				F.clear();
				continue;
			}
            if (beta_v == beta_ && beta_sub > beta_v && alpha_whether_beta_add >= current_alpha) {
                beta_ += 1;
				g.left_index[u][current_alpha] = beta_;
				g.right_index[v][beta_] = current_alpha;
                C.push_back(current_alpha);
                C.push_back(beta_);
                g.is_node_change_in_pair_v1[u][C] = 1;
                C.clear();
                //have_update_u[current_alpha] = beta_;
            }
            else{
                if (skyline_index_u[u]->vector_set[i_iter + 1] < beta_){
                    C.push_back(current_alpha);
                    C.push_back(beta_);
					g.left_index[u][current_alpha] = beta_;
                    //skyline_index_u[u]->vector_set[i_iter + 1] = beta_;
                    g.is_node_change_in_pair_v1[u][C] = 1;
                    C.clear();
                }
				visited_nodes_left.clear(); visited_nodes_right.clear();
				candidate_nodes_left.clear(); candidate_nodes_right.clear();
				visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
				dyn_update_insertion_dfs(g, u, current_alpha, beta_, true, visited_nodes_left, visited_nodes_right,
					candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
				int alpha = current_alpha;
				while(g.left_index[u][alpha] < beta_ + 1 && alpha > 0){
					for (auto i : candidate_nodes_left){
						if(g.left_index[i][alpha] < beta_ + 1) g.left_index[i][alpha] = beta_ + 1;
					}
					for (auto i : candidate_nodes_right){
						if(g.right_index[i][beta_ + 1] < alpha) g.right_index[i][beta_ + 1] = alpha;
					}
					alpha--;
				}
            }
			if(g.left_index[u][current_alpha] == beta_ + 1){
				have_update[F] = 1;
				F.clear();
			}
        }

		for (auto nbr_node : nbr_v){
			B.push_back(compute_ualpha_regarding_specific_beta(g, nbr_node, g.degree_v2[v]));
		}
		int alpha_2 = find_kth(B, B.size(), g.degree_v2[v]);
		B.clear();
		visited_nodes_left.clear(); visited_nodes_right.clear();
		candidate_nodes_left.clear(); candidate_nodes_right.clear();
		visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
		g.right_index[v][g.degree_v2[v]] = max(alpha_2,1);
		dyn_update_insertion_dfs(g, v, alpha_2 + 1, g.degree_v2[v] - 1, false, visited_nodes_left, visited_nodes_right,
			candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		//if(candidate_nodes_right.size() != 1){
			for (auto i : candidate_nodes_left){
				if(g.left_index[i][alpha_2 + 1] < g.degree_v2[v]) g.left_index[i][alpha_2 + 1] = g.degree_v2[v];
			}
			for (auto i : candidate_nodes_right){
				if(g.right_index[i][g.degree_v2[v]] < alpha_2 + 1)g.right_index[i][g.degree_v2[v]] = alpha_2 + 1;
			}
		//}
        B.clear();

		for(int i_iter = 1; i_iter <= skyline_index_v[v]->vector_set.size() - 1; i_iter+=2){
            int current_beta = skyline_index_v[v]->vector_set[i_iter];
            int current_alpha = skyline_index_v[v]->vector_set[i_iter-1];
			if( g.left_index[u].size() < current_alpha) break;
			if (current_alpha > g.degree_v1[u]) continue;
			if (current_beta >= g.left_index[u][current_alpha]){
				if (g.right_index[v][current_beta+1] <= compute_ualpha_regarding_specific_beta(g,u,current_beta+1) && skyline_index_v[v]->vector_set[i_iter - 2] != current_beta + 1){
					for (auto nbr_node : nbr_v){
						B.push_back(compute_ualpha_regarding_specific_beta(g,nbr_node,current_beta+1));
					}
					int alpha = find_kth(B, B.size(),current_beta + 1);
					if(alpha >= skyline_index_v[v]->vector_set[i_iter-3]){
						if (g.right_index[v][current_beta+1] < alpha) {
							g.right_index[v][current_beta+1] = alpha;
							visited_nodes_left.clear(); visited_nodes_right.clear();
							candidate_nodes_left.clear(); candidate_nodes_right.clear();
							visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
							dyn_update_insertion_dfs(g, v, alpha + 1, current_beta, false, visited_nodes_left, visited_nodes_right,
								candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
							for (auto i : candidate_nodes_left){
								g.left_index[i][alpha + 1] = current_beta + 1;
							}
							for (auto i : candidate_nodes_right){
								g.right_index[i][current_beta + 1] = alpha + 1;
							}
						}
					}
					B.clear();
				}
				continue;
			}
            vector<int> current_alpha_and_beta_plus;
            current_alpha_and_beta_plus.push_back(current_alpha);
            current_alpha_and_beta_plus.push_back(current_beta + 1);
			int value;
            for (auto nbr_node : nbr_v){
				value = compute_ualpha_regarding_specific_beta(g, nbr_node,current_beta + 1);
                B.push_back(compute_ualpha_regarding_specific_beta(g, nbr_node, current_beta));
				D.push_back(value);
				if (nbr_node != u) E.push_back(value);
            }
			int beta_v_whether_beta_add_ = find_kth(E, E.size(),current_beta);
            E.clear();
			int alpha_3 = find_kth(D, D.size(),current_beta + 1);
			D.clear();
            int alpha_ = find_kth(B, B.size(),current_beta);
            if (alpha_ == 0) alpha_+= 1;
			int alpha_u = compute_ualpha_regarding_specific_beta(g, u, current_beta);
            int alpha_sub;
            if (current_beta != 1) alpha_sub = find_kth(B, B.size(),current_beta - 1);
            else alpha_sub = alpha_u + 1;
            B.clear();
            vector<int> current_alpha_plus_and_beta;
            current_alpha_plus_and_beta.push_back(alpha_ + 1);
            current_alpha_plus_and_beta.push_back(current_beta);

            if (g.is_node_change_in_pair_v2[v][current_alpha_plus_and_beta] && !g.is_node_change_in_pair_v2[v][current_alpha_and_beta_plus])continue;
            for (auto nbr_node : nbr_u){
                //B.push_back(find_skyline(skyline_index_v, nbr_node,alpha_ + 1, 1, g.num_v2, g.store_skyline_index_v));
				B.push_back(compute_vbeta_regarding_specific_alpha(g, nbr_node, alpha_ + 1));
				if (nbr_node != v) D.push_back(g.right_index[nbr_node][current_beta + 1]);
            }
            int beta_whether_alpha_add = find_kth(B, B.size(),alpha_);
            B.clear();
            int beta_v_whether_beta_add;
            if (D.size() !=0 || current_alpha - 1 != 0) beta_v_whether_beta_add = find_kth(D, D.size(),current_alpha - 1); 
            else beta_v_whether_beta_add = current_alpha;
            D.clear();
            // v-beta change
			F.push_back(alpha_3 + 1);
			F.push_back(current_beta + 1);
            if (!g.is_node_change_in_pair_v2[v][current_alpha_and_beta_plus] && have_update.count(F) == 0){
                //if (current_beta != g.degree_v2[v]){
                    //if(alpha_u >= current_alpha && find_skyline(skyline_index_u, u,current_beta + 1, 0, g.num_v1, g.store_skyline_index_u_beta) >= current_alpha - 1){
					if(alpha_u >= alpha_3 && compute_ualpha_regarding_specific_beta(g, u, current_beta + 1) >= alpha_3 - 1){
						if (beta_v_whether_beta_add >= alpha_3 + 1 && beta_v_whether_beta_add_ >= alpha_3 + 1){ //�Ż� alpha_u + 1 >= current_alpha -> alpha_u > current_alpha
							C.push_back(alpha_3 + 1);
                            C.push_back(current_beta + 1);
							g.right_index[v][current_beta + 1] = alpha_3 + 1;
                           // skyline_index_v[v]->vector_set[i_iter] = current_beta + 1;
                            g.is_node_change_in_pair_v2[v][C] = 1;
                            //have_update_v[current_alpha] = current_beta + 1;
                            C.clear();
                        }
                        else{
							g.right_index[v][current_beta + 1] = alpha_3;
							visited_nodes_left.clear(); visited_nodes_right.clear();
							candidate_nodes_left.clear(); candidate_nodes_right.clear();
							visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
							dyn_update_insertion_dfs(g, v, alpha_3 + 1, current_beta, false, visited_nodes_left, visited_nodes_right,
								candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
							for (auto i : candidate_nodes_left){
								int alpha = alpha_3 + 1;
								while(alpha > 0 && g.left_index[i][alpha] < current_beta + 1){
									g.left_index[i][alpha] = current_beta + 1;
									alpha--;
								}
							}
							for (auto i : candidate_nodes_right){
								int beta = current_beta + 1;
								while(beta > 0 && g.right_index[i][beta] < alpha_3 + 1){
									g.right_index[i][beta] = alpha_3 + 1;
									beta--;
								}
							}
                        }
                    }
                //}
            }
			B.clear();
			if(g.right_index[v][current_beta + 1] < alpha_3) g.right_index[v][current_beta + 1] = alpha_3;
			if(g.right_index[v][current_beta + 1] == alpha_3 + 1){
				have_update[F] = 1;
				F.clear();
			}
            // v-alpha change
			if(current_beta != g.degree_v2[v] - 1){
				C.push_back(alpha_);
				C.push_back(current_beta);
				F.push_back(alpha_ + 1);
				F.push_back(current_beta);
				if (!g.is_node_change_in_pair_v2[v][current_alpha_plus_and_beta] && have_update.count(F) == 0){
					if (alpha_u == alpha_ && alpha_sub > alpha_u && beta_whether_alpha_add >= current_beta) {
						alpha_ += 1;
						//C.push_back(alpha_);
						//C.push_back(current_beta);
						g.right_index[v][current_beta] = alpha_;
						g.left_index[u][alpha_] = current_beta;
						//skyline_index_v[v]->vector_set[i_iter - 1] = alpha_;
						g.is_node_change_in_pair_v2[v][C] = 1;
						//have_update_v[alpha_] = current_beta;
						C.clear();
					}
					else{
						//if (skyline_index_v[v]->vector_set[i_iter - 1] < alpha_){
						if (g.right_index[v][current_beta] < alpha_){
							g.right_index[v][current_beta] = alpha_;
							//skyline_index_v[v]->vector_set[i_iter - 1] = alpha_;
							g.is_node_change_in_pair_v2[v][C] = 1;
							//have_update_v[alpha_] = current_beta;
						}
						visited_nodes_left.clear(); visited_nodes_right.clear();
						candidate_nodes_left.clear(); candidate_nodes_right.clear();
						visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
						dyn_update_insertion_dfs(g, v, alpha_ + 1, current_beta - 1, false, visited_nodes_left, visited_nodes_right,
							candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
						int beta = current_beta;
						while(g.right_index[v][beta] < alpha_ + 1 && beta > 0){
							for (auto i : candidate_nodes_left){
								if(g.left_index[i][alpha_ + 1] < beta) g.left_index[i][alpha_ + 1] = beta;
							}
							for (auto i : candidate_nodes_right){
								if(g.right_index[i][beta] < alpha_ + 1) g.right_index[i][beta] = alpha_ + 1;
							}
							beta--;
						}
					}
					C.clear();
				}
				if(g.right_index[v][current_beta] == alpha_ + 1){
					have_update[F] = 1;
					F.clear();
				}
			}
            current_alpha_and_beta_plus.clear();
            current_alpha_plus_and_beta.clear();
        }
        auto end = chrono::system_clock::now();
		//build_bicore_index(g, bicore_index_u, bicore_index_v);
        //build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}
}