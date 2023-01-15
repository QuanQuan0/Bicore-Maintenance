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

	void dyn_update_removal_dfs(BiGraph& g, vid_t node, int alpha, int beta, bool left_side,
		std::unordered_set<vid_t>& visited_nodes_left, std::unordered_set<vid_t>& visited_nodes_right,
		std::unordered_set<vid_t>& candidate_nodes_left, std::unordered_set<vid_t>& candidate_nodes_right,
		std::unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_left, std::unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_right);

    double update_skyline_index_swap_with_bfs_delete(BiGraph& g, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v,
                                                     std::vector<skyline_block*>& skyline_index_u_reverse, std::vector<skyline_block*>& skyline_index_v_reverse,vid_t u, vid_t v);

	bool dynamic_ergodic_change_by_satisfy_nbr_number(BiGraph& g, vid_t w, int w_is_v1, int current_alpha, int beta_, std::vector<skyline_block*>& skyline_index_u, std::vector<skyline_block*>& skyline_index_v,std::map<std::vector<int>,int>& candidate_nodes);

	void dynamic_ergodic_change_by_satisfy_nbr_number_delete(BiGraph& g, int alpha, int beta,
        vid_t node, bool left_side, std::unordered_set<vid_t>& visited_nodes_left, std::unordered_set<vid_t>& visited_nodes_right,
		std::unordered_set<vid_t>& candidate_nodes_left, std::unordered_set<vid_t>& candidate_nodes_right,
		std::unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_left, std::unordered_map<vid_t, traversal_info_removal>& visited_nodes_2_info_right);

	double update_bicore_index_batch(BiGraph& g, std::vector<std::vector<bicore_index_block*>>& bicore_index_u, std::vector<std::vector<bicore_index_block*>>& bicore_index_v,
		batch_set I, batch_set R);
}
#endif // !DYN_REBUILD_H
