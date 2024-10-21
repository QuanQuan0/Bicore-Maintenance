#pragma once
#ifndef DYN_REBUILD_H
#define DYN_REBUILD_H
#include "bigraph.h"
namespace Dyn_rebuild
{
	struct traversal_info_insertion {
		int current_supports = 0;
		//		int necessary_sup = 0;
		std::vector<vid_t> neighbors_in_candidates;
	};

	struct traversal_info_removal {
		int current_supports = 0;
		//		int necessary_sup = 0;
		std::vector<vid_t> neighbors_not_in_candidates; // not in candidates but in visited nodes
		std::vector<vid_t> nodes_to_be_visited; // in case this node is selected into candidates after other nodes selected into candidates
	};

	// 定义自定义哈希函数，不能完全避免哈希冲突
	struct pair_hash {
		template <class T1, class T2>
		std::size_t operator() (const std::pair<T1, T2>& pair) const {
			return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
		}
	};

	class TopX {
	public:
		TopX(int x) : x(x) {}

		// 插入新元素
		void insert(int num) {
			if (minHeap.size() < x) {
				minHeap.push(num);  // 堆未满时直接插入
			} else if (num > minHeap.top()) {
				minHeap.pop();       // 弹出堆顶（最小值）
				minHeap.push(num);   // 插入新值
			}
		}

		// 获取第x大的元素
		int get_xth_largest() {
			if (minHeap.size() == x) {
				return minHeap.top();  // 堆顶是第x大的元素
			} else {
				return 0;
			}
		}

		// 清空当前的堆
		void clear() {
			while (!minHeap.empty()) {
				minHeap.pop();  // 清空堆
			}
		}

		// 调整x的大小，清空堆
		void resize(int new_x) {
			if (new_x < 1) {
				throw std::invalid_argument("Size must be at least 1.");
			}
			x = new_x;
			clear();  // 清空堆，等待新的插入
		}

		// 删除堆顶元素，获得下一个更大的元素
		void pop_top() {
			if (!minHeap.empty()) {
				minHeap.pop();  // 移除堆顶（第x大的元素）
				x--;
			} else {
				throw std::runtime_error("Heap is empty.");
			}
		}

		// 判断堆是否为空
		bool is_empty() const {
			return minHeap.empty();
		}


	private:
		int x;  // 要维护的第x大的元素
		std::priority_queue<int, std::vector<int>, std::greater<int>> minHeap;  // 使用最小堆
	};


	void dyn_update_removal_dfs(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		std::unordered_set<vid_t>& visited_nodes_left, std::unordered_set<vid_t>& visited_nodes_right,
		std::unordered_set<vid_t>& candidate_nodes_left, std::unordered_set<vid_t>& candidate_nodes_right,
		std::unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_left, std::unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_right);

	double update_bicore_index_without_limit_swap(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition);

	double update_bicore_index_with_limit_swap(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition);

	double update_bicore_index_with_limit_swap_base_opt(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition);

    //double update_skyline_index_swap_with_bfs_delete(BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v,  vid_t u, vid_t v);

	double update_skyline_index_swap_with_bfs(BiGraph& g, vid_t u, vid_t v, int& skyline_pair);
	double update_skyline_index_swap_with_bfs(BiGraph& g, std::vector<skyline_block *>& skyline_index_u, std::vector<skyline_block *>&skyline_index_v, vid_t u, vid_t v);

    int remove_candidate(BiGraph& g, int additional_node, int node, int is_v1, int current_alpha, int beta_, std::map<std::vector<int>,int>& candidate_nodes, std::map<int,int>& is_in_wait_tobe_visit_u, std::map<int,int>& is_in_wait_tobe_visit_v);

    int find_skyline(std::vector<skyline_block*>& skyline_index_uv, vid_t uv, int alpha_or_beta, int is_alpha, int num, std::vector<int>& store_skyline_index_uv);

    int find_kth(std::vector<int> B, int n, int k);

    void dynamic_ergodic_change_by_satisfy_nbr_number_delete(BiGraph& g, int alpha, int beta, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v,std::vector<std::vector<int>>& candidate_nodes,
                                                             std::queue<std::vector<int>>& wait_to_be_visit);

	double update_bicore_index_batch(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		batch_set I, batch_set R);

	double update_bicore_index_with_limit_swap_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition, int threads);

	double update_bicore_index_with_limit_swap_dfs_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		vid_t u, vid_t v, bool addition, int threads);

	double update_bicore_index_batch_parallel(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		batch_set I, batch_set R, int threads);
}
#endif // !DYN_REBUILD_H
