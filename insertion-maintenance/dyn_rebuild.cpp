#include "dyn_rebuild.h"
#include <cmath>
#include <stack>
using namespace std;
// use dynamic_test to test correctness

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

		vector<pair<pair<vid_t,vid_t>,bool>> stack; // pre cur left_side

		stack.reserve(g.num_v1 + g.num_v2); // 提前分配栈的容量
		stack.emplace_back(make_pair(-1, node), left_side);
		bool start_root = true;

		
		while (!stack.empty()) {

			auto [pre_cur, is_left] = stack.back();
			stack.pop_back();

			vid_t pre = pre_cur.first;
			vid_t node = pre_cur.second;

			if(is_left){
				if (visited_nodes_left.count(node)) continue;
				if (!start_root && !candidate_nodes_right.count(pre)) continue;
				start_root = false;

				candidate_nodes_left.insert(node);
				visited_nodes_left.insert(node);
				auto& neigh = g.neighbor_v1[node];

				for (auto v : neigh) {
					if (g.degree_v2[v] > beta) {
						int vbeta = upgrade_vbeta_insertion(g, v, alpha, beta);

						if (vbeta > beta) {
							visited_nodes_2_info_left[node].current_supports++;
						} else if (vbeta == beta) {
							if (visited_nodes_right.count(v) == 0) {
								visited_nodes_2_info_left[node].current_supports++;
								stack.push_back(make_pair(make_pair(node,v),false));
							} else if (candidate_nodes_right.count(v)) {
								visited_nodes_2_info_left[node].current_supports++;
								visited_nodes_2_info_left[node].neighbors_in_candidates.push_back(v);
								visited_nodes_2_info_right[v].neighbors_in_candidates.push_back(node);
							}
						}
					}
				}
				//
				if (visited_nodes_2_info_left[node].current_supports < alpha) {
					casacade_remove_candidate_nodes(node, alpha, beta, true,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
				}

			}
			else {
				if (visited_nodes_right.count(node)) continue;
				if (!start_root && !candidate_nodes_left.count(pre)) continue;

				start_root = false;
				candidate_nodes_right.insert(node);
				visited_nodes_right.insert(node);
				auto& neigh = g.neighbor_v2[node];

				for (auto u : neigh) {
					if (g.degree_v1[u] >= alpha) {
						int ubeta = g.left_index[u][alpha];
						if (ubeta > beta) {
							visited_nodes_2_info_right[node].current_supports++;
						} else if (ubeta == beta) {
							if (visited_nodes_left.find(u) == visited_nodes_left.end()) {
								visited_nodes_2_info_right[node].current_supports++;
								stack.push_back(make_pair(make_pair(node, u), true));
							} else if (visited_nodes_left.find(u) != visited_nodes_left.end() && candidate_nodes_left.find(u) != candidate_nodes_left.end()) {
								visited_nodes_2_info_right[node].current_supports++;
								visited_nodes_2_info_right[node].neighbors_in_candidates.push_back(u);
								visited_nodes_2_info_left[u].neighbors_in_candidates.push_back(node);
							}
						}
					}
				}

				if (visited_nodes_2_info_right[node].current_supports <= beta) {
					casacade_remove_candidate_nodes(node, alpha, beta, false,
						visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right,
						visited_nodes_2_info_left, visited_nodes_2_info_right);
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
		
		if (left_side) {
			visited_nodes_left.insert(node);
			vector<vid_t>& neigh = g.neighbor_v1[node];

			for (auto v : neigh) {
				if (g.degree_v2[v] >= beta) {

					int vbeta = upgrade_vbeta_removal(g, v, alpha, beta);

					// 使用单次查找缓存结果，避免重复调用find()
                	bool v_not_visited = visited_nodes_right.find(v) == visited_nodes_right.end();
                	bool v_not_candidate = candidate_nodes_right.find(v) == candidate_nodes_right.end();

					if (vbeta > beta) {
						visited_nodes_2_info_left[node].current_supports++;
					}
					else if (vbeta == beta) {
						if (v_not_visited) {
							visited_nodes_2_info_left[node].current_supports++;
							visited_nodes_2_info_left[node].nodes_to_be_visited.push_back(v);
						}
						else if (v_not_visited && v_not_candidate) {
							visited_nodes_2_info_left[node].current_supports++;
							visited_nodes_2_info_left[node].neighbors_not_in_candidates.push_back(v);
							visited_nodes_2_info_right[v].neighbors_not_in_candidates.push_back(node);
						}
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
					bool u_not_visited = visited_nodes_left.find(u) == visited_nodes_left.end();
                	bool u_not_candidate = candidate_nodes_left.find(u) == candidate_nodes_left.end();

					if (ubeta > beta) {
						visited_nodes_2_info_right[node].current_supports++;
					}
					else if (ubeta == beta) {
						if (u_not_visited) {
							visited_nodes_2_info_right[node].current_supports++;
							visited_nodes_2_info_right[node].nodes_to_be_visited.push_back(u);
						}
						// can be simplified
						else if (!u_not_visited && u_not_candidate) {
							visited_nodes_2_info_right[node].current_supports++;
							// this node is doomed! must be the end point of insertion edge otherwise it can not be visited!
							if (g.right_index[node][beta] < alpha || g.degree_v2[node] < beta) continue;
							visited_nodes_2_info_right[node].neighbors_not_in_candidates.push_back(u);
							visited_nodes_2_info_left[u].neighbors_not_in_candidates.push_back(node);
						}
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

int find_kth(vector<int> B, int m, int k) {
	int n = B.size();
	if (k <= 0 || k > n) return 0; // handle invalid k values
	std::nth_element(B.begin(), B.begin() + (n - k), B.end());
	return B[n - k];
}

	void rebuild_skyline_index(BiGraph& g, vector<skyline_block*>& skyline_index_u, vector<skyline_block*>& skyline_index_v){
		skyline_index_u.clear(); skyline_index_v.clear(); 
		vector<vector<bicore_index_block*>> bicore_index_u; vector<vector<bicore_index_block*>> bicore_index_v;
        build_bicore_index(g, bicore_index_u, bicore_index_v);
        build_skyline_index(g, bicore_index_u, bicore_index_v, skyline_index_u, skyline_index_v);

	}

	void dyn_crossUpdate_addition_right(BiGraph& g, int alpha, int k_x, vid_t v) {
		int oldalpha;
		for (int beta = k_x; beta > 0; beta--) {
			oldalpha = g.right_index[v][beta];
			if (oldalpha < alpha) {
				g.right_index[v][beta] = alpha;
			}else {
				break;
			}
		}
	}

	void dyn_crossUpdate_addition_left(BiGraph& g, int beta, int k_x, vid_t u) {
		int oldbeta;
		for (int alpha = k_x; alpha > 0; alpha--) {
			oldbeta = g.left_index[u][alpha];
			if (oldbeta < beta){
				g.left_index[u][alpha] = beta;
			}else{
				break;
			}
		}
	}

	vector<int> findDominantPairs_u(const vector<int>& nums, vector<int>& result) {
		int maxY = INT_MIN;  // 当前的最大y值
		result.reserve(nums.size()); // 预留空间，避免频繁的扩容

		// 从右向左遍历vector
		for (int i = nums.size() - 1; i >= 0; --i) {
			// 如果当前元素比maxY大，说明该元素不被支配
			if (nums[i] > maxY) {
				// 将该元素的值和下标加入结果
				result.insert(result.begin(), {i, nums[i]});
				maxY = nums[i];             // 更新最大y值
			}
		}
		result.shrink_to_fit(); // 将预留空间缩减为实际使用的大小
		return result;
	}

	vector<int> findDominantPairs_v(const vector<int>& nums, vector<int>& result) {
		int maxY = INT_MIN;
		result.reserve(nums.size()); 

		for (int i = nums.size() - 1; i >= 0; --i) {
			if (nums[i] > maxY) {
				result.push_back(nums[i]);
				result.push_back(i);
				maxY = nums[i]; 
			}
		}
		result.shrink_to_fit();
		return result;
	}

	int compute_lower_bound_left(BiGraph& g, vid_t u, int alpha, TopX& topX){
		// 计算left side节点在alpha下beta的lower bound
		topX.resize(alpha);
		vector<vid_t>& neigh = g.neighbor_v1[u];
		int ubeta = g.left_index[u][alpha];
		for (vid_t nbr_v : neigh) {
			// 筛选掉不可能属于topX的nbr
			if (ubeta == 0 || g.right_index[nbr_v][ubeta] >= alpha){
				topX.insert(compute_vbeta_regarding_specific_alpha(g, nbr_v, alpha));
			}
		}
		int bound = topX.get_xth_largest();
		if (ubeta < bound) g.left_index[u][alpha] = bound;
		// cout << "alpha: " << alpha << " bound: " << bound << endl;
		return bound;
	}

	int compute_lower_bound_right(BiGraph& g, vid_t v, int beta, TopX& topX){
		// 计算right side节点在beta下alpha的lower bound
		topX.resize(beta);
		vector<vid_t>& neigh = g.neighbor_v2[v];
		int valpha = g.right_index[v][beta];
		for (vid_t nbr_u : neigh) {
			// 筛选掉不可能属于topX的nbr
			if (valpha == 0 || g.left_index[nbr_u][valpha] >= beta){
				topX.insert(compute_ualpha_regarding_specific_beta(g, nbr_u, beta));
			}
		}
		int bound = topX.get_xth_largest();
		if (valpha < bound) g.right_index[v][beta] = bound;
		return bound;
	}

	void case_alpha_increase_left(BiGraph& g, vid_t u, vid_t v, int current_alpha, int current_beta, unordered_set<std::pair<int, int>, pair_hash>& have_update, unordered_set<int>& have_update_alpha,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
        unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right, 
        unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right, TopX& topX, int& bfs_times){
		int bound = compute_lower_bound_left(g, u, current_alpha + 1, topX);
		// 排除重复更新
		if (have_update.find({current_alpha + 1, bound + 1}) != have_update.end()) return;

		//清空容器
		visited_nodes_left.clear(); visited_nodes_right.clear();
		candidate_nodes_left.clear(); candidate_nodes_right.clear();
		visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
		// 调用 DFS 计算候选集
		bfs_times++;
        dyn_update_insertion_dfs(g, u, current_alpha + 1, bound, true, visited_nodes_left, visited_nodes_right,
            candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);

		// 更新候选集
		for (auto u : candidate_nodes_left) {
			if (g.left_index[u][current_alpha + 1] < bound + 1) g.left_index[u][current_alpha + 1] = bound + 1;
			// dyn_crossUpdate_addition_left(g, bound + 1, current_alpha + 1, u);
		}
		for (auto v : candidate_nodes_right) {
			dyn_crossUpdate_addition_right(g, current_alpha + 1, bound + 1, v);
		}

		have_update.insert({current_alpha + 1, bound + 1});
		have_update_alpha.insert(current_alpha + 1);
	}

	void case_beta_increase_right(BiGraph& g, vid_t u, vid_t v, int current_alpha, int current_beta, unordered_set<std::pair<int, int>, pair_hash>& have_update, unordered_set<int>& have_update_beta,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
        unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right, 
        unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right, TopX& topX, int& bfs_times){
		int bound = compute_lower_bound_right(g, v, current_beta + 1, topX);
		// 排除重复更新
		if (have_update.find({bound + 1, current_beta + 1}) != have_update.end()) return;

		//清空容器
		visited_nodes_left.clear(); visited_nodes_right.clear();
		candidate_nodes_left.clear(); candidate_nodes_right.clear();
		visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
		// 调用 DFS 计算候选集
		bfs_times++;
        dyn_update_insertion_dfs(g, v, bound + 1, current_beta, false, visited_nodes_left, visited_nodes_right,
            candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);

		// 更新候选集
		for (auto v : candidate_nodes_right){
			if (g.right_index[v][current_beta + 1] < bound + 1) g.right_index[v][current_beta + 1] = bound + 1;
			// dyn_crossUpdate_addition_right(g, bound + 1, current_beta + 1, v);
		}
		for (auto u : candidate_nodes_left){
			dyn_crossUpdate_addition_left(g, current_beta + 1, bound + 1, u);
		}

		have_update.insert({{bound + 1, current_beta + 1}});
		have_update_beta.insert(current_beta + 1);
	}

	void case_beta_increase_left(BiGraph& g, vid_t u, vid_t v, int current_alpha, int current_beta, unordered_set<std::pair<int, int>, pair_hash>& have_update, unordered_set<int>& have_update_alpha,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
        unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right, 
        unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right, TopX& topX, int& bfs_times){
		int bound = compute_lower_bound_left(g, u, current_alpha, topX);
		// 排除重复更新
		if (have_update.find({current_alpha, bound + 1}) != have_update.end()) return;

		//清空容器
		visited_nodes_left.clear(); visited_nodes_right.clear();
		candidate_nodes_left.clear(); candidate_nodes_right.clear();
		visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
		// 调用 DFS 计算候选集
		bfs_times++;
        dyn_update_insertion_dfs(g, u, current_alpha, bound, true, visited_nodes_left, visited_nodes_right,
            candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		
		// 更新候选集
		for (auto u : candidate_nodes_left) {
			// if (g.left_index[u][current_alpha] < bound + 1) g.left_index[u][current_alpha] = bound + 1;
			dyn_crossUpdate_addition_left(g, bound + 1, current_alpha, u);
		}
		for (auto v : candidate_nodes_right) {
			dyn_crossUpdate_addition_right(g, current_alpha, bound + 1, v);
		}

		have_update.insert({current_alpha, bound + 1});
		have_update_alpha.insert(current_alpha);
	}

	void case_alpha_increase_right(BiGraph& g, vid_t u, vid_t v, int current_alpha, int current_beta, unordered_set<std::pair<int, int>, pair_hash>& have_update, unordered_set<int>& have_update_beta,
		unordered_set<vid_t>& visited_nodes_left, unordered_set<vid_t>& visited_nodes_right,
        unordered_set<vid_t>& candidate_nodes_left, unordered_set<vid_t>& candidate_nodes_right, 
        unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_left, unordered_map<vid_t, traversal_info_insertion>& visited_nodes_2_info_right, TopX& topX, int& bfs_times){
		int bound = compute_lower_bound_right(g, v, current_beta, topX);
		// 排除重复更新
		if (have_update.find({bound + 1, current_beta}) != have_update.end()) return;

		//清空容器
		visited_nodes_left.clear(); visited_nodes_right.clear();
		candidate_nodes_left.clear(); candidate_nodes_right.clear();
		visited_nodes_2_info_left.clear(); visited_nodes_2_info_right.clear();
		// 调用 DFS 计算候选集
		bfs_times++;
        dyn_update_insertion_dfs(g, v, bound + 1, current_beta - 1, false, visited_nodes_left, visited_nodes_right,
            candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right);
		
		// 更新候选集
		for (auto v : candidate_nodes_right){
			// if (g.right_index[v][current_beta] < bound + 1) g.right_index[v][current_beta] = bound + 1;
			dyn_crossUpdate_addition_right(g, bound + 1, current_beta, v);
		}
		for (auto u : candidate_nodes_left){
			dyn_crossUpdate_addition_left(g, current_beta, bound + 1, u);
		}

		have_update.insert({{bound + 1, current_beta}});
		have_update_beta.insert(current_beta);
	}

	double update_skyline_index_swap_with_bfs(BiGraph& g, vector<skyline_block *>& skyline_index_u, vector<skyline_block *>&skyline_index_v, vid_t u, vid_t v){
		return 0.00;
	}

	double update_skyline_index_swap_with_bfs(BiGraph& g, vid_t u, vid_t v, int& skyline_pair){
		auto start = chrono::system_clock::now();

		int bfs_times = 0;

		g.addEdge(u, v);

		// 生成接下来分析所需的skyline pair
		vector<int> skyline_u, skyline_v;
		findDominantPairs_u(g.left_index[u], skyline_u);
		findDominantPairs_v(g.right_index[v], skyline_v);

		// 初始化
		vector<vid_t>& nbr_u = g.neighbor_v1[u]; 
       		vector<vid_t>& nbr_v = g.neighbor_v2[v];
		int u_deg = g.degree_v1[u];
		int v_deg = g.degree_v2[v];
		// cout << "deg u: " << u_deg << " " << skyline_u.size() << endl;

		//case, u or v's degree = 1
		if(u_deg == 1){
			g.left_index[u].push_back(v_deg);
			g.right_index[v].push_back(u_deg);
			for(auto nbr_node : nbr_v){
				if(g.left_index[nbr_node][1] < v_deg) 
                    g.left_index[nbr_node][1] = v_deg;
			}
			auto end = chrono::system_clock::now();
            return chrono::duration<double>(end - start).count();
		}

		if(v_deg == 1){
			g.right_index[v].push_back(u_deg);
			g.left_index[u].push_back(v_deg);
			for (auto nbr_node : nbr_u){
				if(g.right_index[nbr_node][1] < u_deg) g.right_index[nbr_node][1] = u_deg;
			}
			auto end = chrono::system_clock::now();
            return chrono::duration<double>(end - start).count();
		}

		g.left_index[u].push_back(0);
		g.right_index[v].push_back(0);

		// (0,0)也是一个skyline
		if (skyline_u[0] != 1) skyline_u.insert(skyline_u.begin(), {1, g.left_index[u][1]});
		if (skyline_v.empty() || skyline_v.back() != 1) skyline_v.insert(skyline_v.end(), {g.right_index[v][1], 1});

		unordered_set<std::pair<int, int>, pair_hash> have_update; // 记录防止重复更新
		unordered_set<int> have_update_alpha;
		unordered_set<int> have_update_beta;
		unordered_set<vid_t> visited_nodes_left, visited_nodes_right;
        unordered_set<vid_t> candidate_nodes_left, candidate_nodes_right;
        unordered_map<vid_t, traversal_info_insertion> visited_nodes_2_info_left, visited_nodes_2_info_right;

		// 优先队列，用于找第X大
		TopX topX(1); 

		for(int i_iter = skyline_u.size()-2; i_iter >= 0; i_iter-=2){
			int current_alpha = skyline_u[i_iter];
            int current_beta = skyline_u[i_iter+1];
			// cout << "u: " << u << " " << current_alpha << " " << current_beta << endl;

			// if (current_alpha > current_beta) continue;

			// case: alpha increase
			int bv_alpha_p_one = compute_vbeta_regarding_specific_alpha(g, v, current_alpha + 1);
			int bu_alpha_p_one = g.left_index[u][current_alpha + 1];
			// 更新条件：bv_alpha_p_one >= bu_alpha_p_one；current_alpha+1未被计算过
			if (bv_alpha_p_one >= bu_alpha_p_one && have_update_alpha.find(current_alpha + 1) == have_update_alpha.end()) 
				case_alpha_increase_left(g, u, v, current_alpha, current_beta, have_update, have_update_alpha,
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right, topX, bfs_times);

			if (v_deg < current_beta) break;

			// case: beta increase
			int bv_alpha = compute_vbeta_regarding_specific_alpha(g, v, current_alpha);
			// 更新条件: bv_alpha >= current_beta; current_alpha未被计算过
			if (bv_alpha >= current_beta && have_update_alpha.find(current_alpha) == have_update_alpha.end())
				case_beta_increase_left(g, u, v, current_alpha, current_beta, have_update, have_update_alpha,
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right, topX, bfs_times);
		}

		for (int i_iter = 1; i_iter <= skyline_v.size() - 1; i_iter+=2){
			int current_beta = skyline_v[i_iter];
            int current_alpha = skyline_v[i_iter-1];

			// cout << "v: " << v << " " << current_alpha << " " << current_beta << " " << v_deg << endl;

			// if (current_beta > current_alpha) continue;

			// case: beta increase
			int au_beta_p_one = compute_ualpha_regarding_specific_beta(g, u, current_beta + 1);
			int av_beta_p_one = g.right_index[v][current_beta + 1];
			// 更新条件: au_beta_p_one >= av_beta_p_one; current_beta+1未被计算过
			if (au_beta_p_one >= av_beta_p_one && have_update_beta.find(current_beta + 1) == have_update_beta.end()){
				case_beta_increase_right(g, u, v, current_alpha, current_beta, have_update, have_update_beta,
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right, topX, bfs_times);}
			
			if (u_deg < current_alpha) break;

			// case: alpha increase
			int au_alpha = compute_ualpha_regarding_specific_beta(g, u, current_beta);
			// 更新条件: au_alpha >= current_alpha; current_beta未被计算过
			if (au_alpha >= current_alpha && have_update_beta.find(current_beta) == have_update_beta.end())
				case_alpha_increase_right(g, u, v, current_alpha, current_beta, have_update, have_update_beta, 
					visited_nodes_left, visited_nodes_right, candidate_nodes_left, candidate_nodes_right, visited_nodes_2_info_left, visited_nodes_2_info_right, topX, bfs_times);
		}
		cout << "bfs times: " << bfs_times << endl;
		auto end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		return elapsed_seconds.count();
	}
}
